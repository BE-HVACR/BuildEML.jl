# Performance data of the DX coils: curves, stage records and the capacity block.
#
# Refer to the MBL v11:
# - biquadratic: Buildings.Utilities.Math.Functions.biquadratic
# - extended_polynomial: Buildings.Fluid.Utilities.extendedPolynomial
# - DX_CURVE_I, DX_CURVE_II, DX_CURVE_III: Buildings.Fluid.DXSystems.Cooling.AirSource.Examples.PerformanceCurves
#   Curve_I, Curve_II, Curve_III (Curve_II fitted from the EnergyPlus 7.1 example file
#   DXCoilSystemAuto.idf)
# - dx_stage: Buildings.Fluid.DXSystems.Cooling.AirSource.Data.Generic.BaseClasses Stage
#   and NominalValues
# - DXCapacityAirSource: Buildings.Fluid.DXSystems.BaseClasses.CapacityAirSource with nSta = 1

# x1 and x2 in degC
biquadratic(a, x1, x2) = a[1] + x1 * (a[2] + a[3] * x1) + x2 * (a[4] + a[5] * x2) + a[6] * x1 * x2

# the polynomial sum(c[i] * x^(i-1)) inside [xmin, xmax], extended linearly
# outside with a continuous first derivative
function extended_polynomial(x, c, xmin, xmax)
    poly(z) = sum(c[i] * z^(i - 1) for i in eachindex(c))
    dpoly(z) = sum((i - 1) * c[i] * z^(i - 2) for i in 2:lastindex(c))
    return ifelse(x < xmin, poly(xmin) + (x - xmin) * dpoly(xmin),
        ifelse(x > xmax, poly(xmax) + (x - xmax) * dpoly(xmax),
            poly(x)))
end

const DX_CURVE_I = (
    capFunT = [0.942587793, 0.009543347, 0.000683770, -0.011042676, 0.000005249, -0.000009720],
    capFunFF = [0.8, 0.2, 0.0, 0.0],
    EIRFunT = [0.342414409, 0.034885008, -0.000623700, 0.004977216, 0.000437951, -0.000728028],
    EIRFunFF = [1.1552, -0.1808, 0.0256, 0.0],
    TConInMin = 291.15, TConInMax = 319.26111,
    TEvaInMin = 285.92778, TEvaInMax = 297.03889,
    ffMin = 0.6, ffMax = 1.8,
)

const DX_CURVE_II = (
    capFunT = [0.766956, 0.0107756, -0.0000414703, 0.00134961, -0.000261144, 0.000457488],
    capFunFF = [0.8, 0.2, 0.0, 0.0],
    EIRFunT = [0.297145, 0.0430933, -0.000748766, 0.00597727, 0.000482112, -0.000956448],
    EIRFunFF = [1.1552, -0.1808, 0.0256, 0.0],
    TConInMin = 291.15, TConInMax = 319.26111,
    TEvaInMin = 285.92778, TEvaInMax = 297.03889,
    ffMin = 0.5, ffMax = 1.5,
)

const DX_CURVE_III = (
    capFunT = [0.476428, 0.0401147, 0.000226411, -0.000827136, -0.0000073224, -0.000446278],
    capFunFF = [0.47278589, 1.2433415, -1.0387055, 0.32257813],
    EIRFunT = [0.632475, -0.0121321, 0.000507773, 0.0155377, 0.00027284, -0.000679201],
    EIRFunFF = [1.0079484, 0.34544129, -0.6922891, 0.33889943],
    TConInMin = 297.03889, TConInMax = 319.26111,
    TEvaInMin = 285.92778, TEvaInMax = 297.03889,
    ffMin = 0.5, ffMax = 1.5,
)

"""
    dx_stage(; Qflow_nominal, COP_nominal, SHR_nominal, mflow_nominal,
               perf_curves = DX_CURVE_II, spe = 1.0, TEvaIn_nominal = 292.55,
               phiIn_nominal = 0.5, p_nominal = 101325.0)

Performance data of one DX coil stage, as in MBL's `Stage` and `NominalValues`
records. `spe` is the compressor speed of the stage [1/s]; it only matters for
variable-speed coils, where the stages are interpolated by speed.
"""
dx_stage(; Qflow_nominal, COP_nominal, SHR_nominal, mflow_nominal,
    perf_curves = DX_CURVE_II, spe = 1.0, TEvaIn_nominal = 292.55,
    phiIn_nominal = 0.5, p_nominal = 101325.0) =
    (; Qflow_nominal, COP_nominal, SHR_nominal, mflow_nominal,
        perf_curves, spe, TEvaIn_nominal, phiIn_nominal, p_nominal)

dx_mflow_nominal_last(stages) = stages[end].mflow_nominal
dx_speRat_first(stages) = stages[1].spe / stages[end].spe

"""
Cooling capacity and EIR of a single-stage DX coil from performance curves.

Ports: `on` (0 off, 1 on), `TConIn` (condenser inlet [K]), `TEvaIn`
(evaporator air inlet [K], wet bulb for a wet coil, dry bulb for a dry coil),
`mflow` (evaporator air mass flow [kg/s]), `Qflow` (total cooling capacity
[W], <= 0 when on), `EIR` (energy input ratio, 1/COP).

`perf_curves` packs the curve coefficients and their validity ranges;
the default is MBL's `Curve_II` (`DX_CURVE_II`). `Qflow_nominal`,
`COP_nominal` and `mflow_nominal` are the rating-point values the curves
are normalized at.
"""
@mtkmodel DXCapacityAirSource begin
    @structural_parameters begin
        perf_curves = DX_CURVE_II
    end

    @parameters begin
        Qflow_nominal = -21000.0            # [W] < 0
        COP_nominal   = 3.0                 # [-]
        mflow_nominal = 1.5                 # [kg/s]
        mflow_small   = 1e-4 * mflow_nominal   # [kg/s] regularization
    end

    @components begin
        on     = RealInput()
        TConIn = RealInput()    # [K]
        TEvaIn = RealInput()    # [K]
        mflow  = RealInput()    # [kg/s]

        Qflow = RealOutput()    # [W]
        EIR   = RealOutput()    # [-]
    end

    @variables begin
        ff(t)        # [-] flow fraction
        cap_T(t), [guess = 1.0]
        cap_FF(t), [guess = 1.0]
        EIR_T(t), [guess = 1.0]
        EIR_FF(t), [guess = 1.0]
        corFac(t), [guess = 1.0]
        Qflow_(t), [irreducible = true, guess = 0.0]   # [W]
        EIR_(t), [irreducible = true, guess = 0.3]     # [-]
    end

    @equations begin
        ff ~ reg_max(mflow.u, mflow_small, mflow_small / 10) / mflow_nominal

        cap_T ~ reg_max(
            biquadratic(perf_curves.capFunT, TEvaIn.u - 273.15, TConIn.u - 273.15),
            0.001, 0.0001)
        cap_FF ~ extended_polynomial(ff, perf_curves.capFunFF, perf_curves.ffMin, perf_curves.ffMax)

        EIR_T ~ reg_max(
            biquadratic(perf_curves.EIRFunT, TEvaIn.u - 273.15, TConIn.u - 273.15),
            0.001, 0.0001)
        EIR_FF ~ extended_polynomial(ff, perf_curves.EIRFunFF, perf_curves.ffMin, perf_curves.ffMax)

        corFac ~ reg_heaviside(ff - perf_curves.ffMin / 4, max(1e-15, perf_curves.ffMin / 8))

        Qflow_ ~ on.u * corFac * cap_T * cap_FF * Qflow_nominal
        EIR_   ~ on.u * corFac * EIR_T * EIR_FF / COP_nominal
        Qflow.u ~ Qflow_
        EIR.u   ~ EIR_
    end
end
