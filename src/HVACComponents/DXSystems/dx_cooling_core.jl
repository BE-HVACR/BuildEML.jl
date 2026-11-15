# Internal parts of the DX cooling coils: wet/dry coil condition with one or
# more stages.
#
# Refer to the MBL v11 (Buildings.Fluid.DXSystems.Cooling.BaseClasses):
# - dx_uacp: UACp, solved at the rating point when the model is built
# - speed_shift: Functions.speedShift, with Buildings.Utilities.Math.Functions
#   splineDerivatives, isMonotonic and cubicHermiteLinearExtrapolation
# - DXCoolingCore: DXCooling with variableSpeedCoil = true, assembled from WetCoil,
#   DryCoil (DXCapacityAirSource per stage), SpeedShift, ApparatusDewPoint,
#   SensibleHeatRatio, Condensation and DryWetSelector
# The connectors use the library's dry-air basis (w, dry-air mass flow), while
# the intermediate quantities follow MBL's moist-air basis (mass fraction X,
# moist-air mass flow, enthalpy per kg of moist air).
# The accumulated-water evaporation model (Evaporation) is not implemented.

X_of_w(w) = w / (1 + w)
w_of_X(X) = X / (1 - X)
h_moist(T, w) = h_Tw(T, w) / (1 + w)   # J/kg moist air

"""
    dx_uacp(Qflow_nominal, SHR_nominal, mflow_nominal,
            TEvaIn_nominal = 292.55, phiIn_nominal = 0.5, p_nominal = 101325.0)

UA/cp of a DX coil [kg/s], calibrated at the rating point from the nominal
sensible heat ratio, as in MBL's UACp. The apparatus dew point is the first
intersection below the outlet temperature of the inlet-outlet line with the
saturation curve, found by a downward scan and bisection.
"""
function dx_uacp(Qflow_nominal, SHR_nominal, mflow_nominal,
    TEvaIn_nominal = 292.55, phiIn_nominal = 0.5, p_nominal = 101325.0)
    k = 0.621964713077499
    X_pSatpphi(psat, phi) = phi * k / (k * phi + p_nominal / psat - phi)

    T = TEvaIn_nominal
    X_in = X_pSatpphi(p_saturation(T), phiIn_nominal)
    h_in = h_moist(T, w_of_X(X_in))
    h_out = h_in + Qflow_nominal / mflow_nominal
    Cp = cp_da * (1 - X_in) + cp_vapor * X_in
    T_out = T + SHR_nominal * Qflow_nominal / (mflow_nominal * Cp)
    X_out = X_in + (h_out - h_in) * (1 - SHR_nominal) / h_g(T)

    X_sat(Tadp) = X_pSatpphi(p_saturation(Tadp), 1.0)
    res(Tadp) = (X_in - X_out) * (T - Tadp) - (X_in - X_sat(Tadp)) * (T - T_out)
    hi = T_out
    lo = hi - 0.5
    while sign(res(lo)) == sign(res(hi))
        hi = lo
        lo -= 0.5
        lo < 223.15 && error("dx_uacp: no apparatus dew point found above -50 degC")
    end
    for _ in 1:100
        mid = (lo + hi) / 2
        sign(res(mid)) == sign(res(hi)) ? (hi = mid) : (lo = mid)
        hi - lo < 1e-10 && break
    end
    T_adp = (lo + hi) / 2

    h_adp = h_moist(T_adp, w_of_X(X_sat(T_adp)))
    bypass = reg_limit((h_out - h_adp) / (h_in - h_adp), 1e-3, 0.999, 1e-4)
    return -mflow_nominal * log(bypass)
end

"""
    speed_shift(spe, speSet, u)

Value at compressor speed `spe` interpolated from the stage values `u` at the
stage speeds `speSet`, with an extra point (0, 0) at zero speed, as in MBL's
speedShift: a cubic Hermite spline whose derivatives are limited to keep a
monotonic sequence monotonic, extrapolated linearly outside the data. Works
on numbers and on symbolic expressions.
"""
function speed_shift(spe, speSet, u)
    x = vcat(0.0, speSet)
    y = vcat(0.0, collect(u))
    n = length(x)

    delta = [(y[i + 1] - y[i]) / (x[i + 1] - x[i]) for i in 1:(n - 1)]
    if n == 2
        d = [delta[1], delta[1]]
    else
        d = vcat(delta[1], [(delta[i - 1] + delta[i]) / 2 for i in 2:(n - 1)], delta[n - 1])
        inc = prod(ifelse(y[i] <= y[i + 1], 1.0, 0.0) for i in 1:(n - 1))
        dec = prod(ifelse(y[i] >= y[i + 1], 1.0, 0.0) for i in 1:(n - 1))
        mono = ifelse(y[1] >= y[n], dec, inc) > 0.5
        for i in 1:(n - 1)
            flat = abs(delta[i]) < 1e-60
            alpha = d[i] / delta[i]
            beta = d[i + 1] / delta[i]
            s = alpha^2 + beta^2
            di = ifelse(flat, 0.0, ifelse(s > 9, delta[i] * alpha * 3 / sqrt(s), d[i]))
            dj = ifelse(flat, 0.0, ifelse(s > 9, delta[i] * beta * 3 / sqrt(s), d[i + 1]))
            d[i] = ifelse(mono, di, d[i])
            d[i + 1] = ifelse(mono, dj, d[i + 1])
        end
    end

    function hermite_lin(i)
        x1, x2, y1, y2, d1, d2 = x[i], x[i + 1], y[i], y[i + 1], d[i], d[i + 1]
        h = x2 - x1
        tt = (spe - x1) / h
        yh = y1 * (2tt^3 - 3tt^2 + 1) + h * d1 * (tt^3 - 2tt^2 + tt) +
             y2 * (-2tt^3 + 3tt^2) + h * d2 * (tt^3 - tt^2)
        return ifelse(spe <= x1, y1 + (spe - x1) * d1,
            ifelse(spe >= x2, y2 + (spe - x2) * d2, yh))
    end
    res = hermite_lin(1)
    for j in 2:(n - 1)
        res = ifelse(spe > x[j], hermite_lin(j), res)
    end
    return res
end

"""
Wet/dry coil condition of a DX cooling coil (signal level, no air ports).

Inputs: `on` (0 off, 1 on), `speRat` (compressor speed ratio [0, 1], speed =
`speRat` times the speed of the last stage), `TConIn` [K], `T_in`/`w_in`/`p_in`
(evaporator air inlet state), `mflow` [kg_da/s]. Outputs: `Qflow` [W] (total,
<= 0), `EIR` [-], `SHR` [-], `mWat_flow` [kg/s] (condensed water, <= 0).

`stages` is a vector of `dx_stage` records. Each stage evaluates its capacity
curves, and the capacity, EIR and UA/cp are interpolated by speed with
`speed_shift`. The wet coil evaluates the curves at the inlet wet-bulb
temperature and resolves the apparatus dew point from UA/cp, the dry coil
evaluates them at the dry-bulb temperature, and the outputs blend smoothly
between the two based on X_adp - X_in. UA/cp of each stage is calibrated from
its nominal values when the model is built.
"""
@mtkmodel DXCoolingCore begin
    @structural_parameters begin
        stages = [dx_stage(Qflow_nominal = -21000.0, COP_nominal = 3.0,
            SHR_nominal = 0.8, mflow_nominal = 1.5)]
        N = length(stages)
        speSet = [s.spe for s in stages]
        uacp = [dx_uacp(s.Qflow_nominal, s.SHR_nominal, s.mflow_nominal,
            s.TEvaIn_nominal, s.phiIn_nominal, s.p_nominal) for s in stages]
    end

    @parameters begin
        mflow_small = 1e-4 * dx_mflow_nominal_last(stages)   # [kg/s]
        h_adp_min = h_moist(276.15, w_of_X(0.004667))     # [J/kg]
    end

    @components begin
        on     = RealInput()
        speRat = RealInput()
        TConIn = RealInput()    # [K]
        T_in   = RealInput()    # [K]
        w_in   = RealInput()    # [kg_w/kg_da]
        p_in   = RealInput()    # [Pa]
        mflow  = RealInput()    # [kg_da/s]

        wetCap = [DXCapacityAirSource(
            perf_curves = stages[i].perf_curves, Qflow_nominal = stages[i].Qflow_nominal,
            COP_nominal = stages[i].COP_nominal, mflow_nominal = stages[i].mflow_nominal)
                  for i in 1:N]
        dryCap = [DXCapacityAirSource(
            perf_curves = stages[i].perf_curves, Qflow_nominal = stages[i].Qflow_nominal,
            COP_nominal = stages[i].COP_nominal, mflow_nominal = stages[i].mflow_nominal)
                  for i in 1:N]

        Qflow     = RealOutput()    # [W]
        EIR       = RealOutput()    # [-]
        SHR       = RealOutput()    # [-]
        mWat_flow = RealOutput()    # [kg/s]
    end

    @variables begin
        spe(t)                          # [1/s] compressor speed
        Q_wet(t), [irreducible = true]
        EIR_wet(t), [irreducible = true]
        Q_dry(t), [irreducible = true]
        EIR_dry(t), [irreducible = true]
        UAcp(t), [guess = 1.0]          # [kg/s]
        T_wb(t), [guess = 288.15]       # [K] inlet wet bulb
        X_in(t), [guess = 0.008]        # [kg_w/kg] inlet mass fraction
        m_moist(t), [guess = 1.0]       # [kg/s] moist-air mass flow
        h_in(t), [guess = 40000.0]      # [J/kg]
        mflow_nz(t), [guess = 1.0]
        bypass(t), [guess = 0.25]
        delta_h(t), [guess = 20000.0]
        h_adp(t), [guess = 30000.0]
        T_adp(t), [guess = 283.15]
        w_adp(t), [guess = 0.007]
        shr_wet(t), [guess = 0.8]
        fraDry(t), [guess = 0.0]
    end

    @equations begin
        0 ~ wetbulb_balance_Twp(T_wb, T_in.u, w_in.u, p_in.u)

        X_in ~ X_of_w(w_in.u)
        m_moist ~ mflow.u * (1 + w_in.u)
        h_in ~ h_moist(T_in.u, w_in.u)
        spe ~ speRat.u * speSet[end]

        [wetCap[i].on.u ~ on.u for i in 1:N]...
        [wetCap[i].TConIn.u ~ TConIn.u for i in 1:N]...
        [wetCap[i].TEvaIn.u ~ T_wb for i in 1:N]...
        [wetCap[i].mflow.u ~ m_moist for i in 1:N]...
        [dryCap[i].on.u ~ on.u for i in 1:N]...
        [dryCap[i].TConIn.u ~ TConIn.u for i in 1:N]...
        [dryCap[i].TEvaIn.u ~ T_in.u for i in 1:N]...
        [dryCap[i].mflow.u ~ m_moist for i in 1:N]...

        Q_wet ~ speed_shift(spe, speSet, [wetCap[i].Qflow_ for i in 1:N])
        EIR_wet ~ speed_shift(spe, speSet, [wetCap[i].EIR_ for i in 1:N])
        Q_dry ~ speed_shift(spe, speSet, [dryCap[i].Qflow_ for i in 1:N])
        EIR_dry ~ speed_shift(spe, speSet, [dryCap[i].EIR_ for i in 1:N])
        UAcp ~ reg_max(speed_shift(spe, speSet, uacp), uacp[end] / 1e3, uacp[end] / 1e4)

        mflow_nz ~ reg_max(m_moist, mflow_small, 0.1 * mflow_small)
        bypass ~ reg_limit(exp(-UAcp / mflow_nz), 0.01, 0.99, 0.001)
        delta_h ~ reg_min(-Q_wet / mflow_nz / (1 - bypass), 0.999 * h_in, 0.0001)
        h_adp ~ reg_min(reg_max(h_in - delta_h, h_adp_min, 10.0), h_in + 100.0, 10.0)

        w_adp ~ 0.62198 * p_saturation(T_adp) / (p_in.u - p_saturation(T_adp))
        h_adp ~ h_moist(T_adp, w_adp)

        shr_wet ~ on.u * reg_min(
            (h_moist(T_in.u, w_adp) - h_adp) / reg_max(100.0, h_in - h_adp, 10.0),
            0.999, 0.0001)

        fraDry ~ splice_fun(1.0, 0.0, (X_of_w(w_adp) - X_in) + 0.0001, 0.0001)

        Qflow.u ~ fraDry * Q_dry + (1 - fraDry) * Q_wet
        EIR.u ~ fraDry * EIR_dry + (1 - fraDry) * EIR_wet
        SHR.u ~ fraDry + (1 - fraDry) * shr_wet
        mWat_flow.u ~ (1 - fraDry) * (1 - shr_wet) * Q_wet / h_fg
    end
end
