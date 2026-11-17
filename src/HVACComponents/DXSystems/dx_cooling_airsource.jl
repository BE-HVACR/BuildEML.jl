# Air-source DX cooling coils.
#
# Refer to the MBL v11 Buildings.Fluid.DXSystems.Cooling.AirSource:
# - DXCoolingSingleSpeed: SingleSpeed
# - DXCoolingVariableSpeed: VariableSpeed
# Both with computeReevaporation = false. Reevaporation of the water accumulated
# on the coil after the compressor stops is not modeled yet.

"""
Single-speed air-source DX cooling coil.

Ports: `port_a` (inlet), `port_b` (outlet), `on` (0 off, 1 on), `TConIn`
(condenser inlet, outdoor dry bulb [K]). Outputs: `Qflow` (total cooling [W],
<= 0), `QSen_flow`, `QLat_flow`, `P` (electric power [W]).

Air passes a pressure drop (`dp_nominal` at the moist-air `mflow_nominal`)
and a mixing volume of `mflow_nominal * tau / rho_air_ref`. `DXCoolingCore`
evaluates the coil from the inlet state, and its total cooling rate and
condensed water are injected into the volume.
"""
@mtkmodel DXCoolingSingleSpeed begin
    @structural_parameters begin
        perf_curves = DX_CURVE_II
        Qflow_nominal = -21000.0
        COP_nominal = 3.0
        mflow_nominal = 1.5
        SHR_nominal = 0.8
        TEvaIn_nominal = 292.55
        phiIn_nominal = 0.5
        p_nominal = 101325.0
        if_steady_state::Bool = false
    end

    @parameters begin
        dp_nominal = 100.0   # [Pa]
        tau = 30.0           # [s]
        deltaM = 0.3         # [-]
    end

    @components begin
        port_a = AirPort()
        port_b = AirPort()
        on     = RealInput()
        TConIn = RealInput()    # [K]

        preDro = AirPressureDrop(
            mflow_nominal = mflow_nominal, dp_nominal = dp_nominal, deltaM = deltaM)
        vol = AirMixingVolumeNodeN(
            N = 2, V = mflow_nominal * tau / rho_air_ref,
            if_moisture_input = true, if_steady_state = if_steady_state)
        core = DXCoolingCore(stages = [dx_stage(
            Qflow_nominal = Qflow_nominal, COP_nominal = COP_nominal,
            SHR_nominal = SHR_nominal, mflow_nominal = mflow_nominal,
            perf_curves = perf_curves, TEvaIn_nominal = TEvaIn_nominal,
            phiIn_nominal = phiIn_nominal, p_nominal = p_nominal)])
        preHea = PrescribedHeatFlow()

        Qflow     = RealOutput()    # [W]
        QSen_flow = RealOutput()    # [W]
        QLat_flow = RealOutput()    # [W]
        P         = RealOutput()    # [W]
    end

    @equations begin
        connect(port_a, preDro.port_a)
        connect(preDro.port_b, vol.airport_1)
        connect(vol.airport_2, port_b)

        connect(on, core.on)
        core.speRat.u ~ 1.0
        connect(TConIn, core.TConIn)
        core.T_in.u ~ preDro.T_in
        core.w_in.u ~ preDro.w_in
        core.p_in.u ~ port_a.p
        core.mflow.u ~ preDro.mflow

        preHea.Q_flow.u ~ core.Qflow.u
        connect(preHea.port, vol.heatport)
        vol.water_mflow.u ~ core.mWat_flow.u

        Qflow.u ~ core.Qflow.u
        QSen_flow.u ~ core.SHR.u * core.Qflow.u
        QLat_flow.u ~ (1 - core.SHR.u) * core.Qflow.u
        P.u ~ -core.EIR.u * core.Qflow.u
    end
end

"""
Variable-speed air-source DX cooling coil.

Ports: `port_a` (inlet), `port_b` (outlet), `speRat` (compressor speed ratio
[0, 1], speed = `speRat` times the speed of the last stage), `TConIn`
(condenser inlet, outdoor dry bulb [K]). Outputs: `Qflow` (total cooling [W],
<= 0), `QSen_flow`, `QLat_flow`, `P` (electric power [W]).

`stages` is a vector of `dx_stage` records ordered by speed; capacity, EIR and
UA/cp are interpolated by speed between them. The pressure drop and the mixing
volume use the nominal air flow of the last stage.

`if_hysteresis_shutoff` sets the behavior at low speed. `speRat_1` below is
the speed ratio of the first stage.

- `false` (default, cycling behavior): below `speRat_1`, the coil cycles at the first stage.
  Total cooling, condensed water and power are the first-stage values times
  the runtime fraction `PLR = speRat / speRat_1`. Within `deltaSpeRat` of
  `speRat_1`, the coil runs at `smooth_max_upper(speRat, speRat_1)` and `PLR`
  is `speRat` divided by it, so `PLR` does not exceed 1. The cycling losses at
  start-up are not modeled yet (EnergyPlus corrects them with a part-load
  fraction curve).
- `true` (MBL hysteresis behavior): the coil runs at the speed `speRat`, and below `speRat_1` the
  performance is interpolated toward zero at zero speed. The coil switches off
  when `speRat` drops below `minSpeRat - speRatDeaBan / 2` and on again above
  `minSpeRat + speRatDeaBan / 2`. This hysteresis relaxes to its target with
  the time constant `tau_on`.

`if_switching_stage_interpolation` is passed to `DXCoolingCore`. With `true`
there and `if_hysteresis_shutoff = true`, the coil reproduces MBL.
"""
@mtkmodel DXCoolingVariableSpeed begin
    @structural_parameters begin
        stages = [dx_stage(Qflow_nominal = -21000.0, COP_nominal = 3.0,
            SHR_nominal = 0.8, mflow_nominal = 1.5)]
        if_hysteresis_shutoff::Bool = false
        if_switching_stage_interpolation::Bool = false
        if_steady_state::Bool = false
    end

    @parameters begin
        dp_nominal = 100.0     # [Pa]
        tau = 30.0             # [s]
        deltaM = 0.3           # [-]
        if if_hysteresis_shutoff
            minSpeRat = 0.2        # [-]
            speRatDeaBan = 0.05    # [-]
            tau_on = 0.01          # [s]
        else
            deltaSpeRat = 0.01     # [-]
        end
    end

    @components begin
        port_a = AirPort()
        port_b = AirPort()
        speRat = RealInput()
        TConIn = RealInput()    # [K]

        preDro = AirPressureDrop(
            mflow_nominal = dx_mflow_nominal_last(stages), dp_nominal = dp_nominal, deltaM = deltaM)
        vol = AirMixingVolumeNodeN(
            N = 2, V = dx_mflow_nominal_last(stages) * tau / rho_air_ref,
            if_moisture_input = true, if_steady_state = if_steady_state)
        core = DXCoolingCore(stages = stages,
            if_switching_stage_interpolation = if_switching_stage_interpolation)
        preHea = PrescribedHeatFlow()

        Qflow     = RealOutput()    # [W]
        QSen_flow = RealOutput()    # [W]
        QLat_flow = RealOutput()    # [W]
        P         = RealOutput()    # [W]
    end

    @variables begin
        if if_hysteresis_shutoff
            on_sta(t), [guess = 0.0]
        else
            PLR(t)    # [-]
        end
    end

    @equations begin
        connect(port_a, preDro.port_a)
        connect(preDro.port_b, vol.airport_1)
        connect(vol.airport_2, port_b)

        connect(TConIn, core.TConIn)
        core.T_in.u ~ preDro.T_in
        core.w_in.u ~ preDro.w_in
        core.p_in.u ~ port_a.p
        core.mflow.u ~ preDro.mflow

        if if_hysteresis_shutoff
            tau_on * D(on_sta) ~ ifelse(speRat.u >= minSpeRat + speRatDeaBan / 2, 1.0,
                ifelse(speRat.u <= minSpeRat - speRatDeaBan / 2, 0.0, on_sta)) - on_sta
            core.on.u ~ on_sta
            core.speRat.u ~ speRat.u
            Qflow.u ~ core.Qflow.u
            vol.water_mflow.u ~ core.mWat_flow.u
        else
            core.on.u ~ 1.0
            core.speRat.u ~ smooth_max_upper(speRat.u, dx_speRat_first(stages); δ = deltaSpeRat)
            PLR ~ speRat.u / core.speRat.u
            Qflow.u ~ PLR * core.Qflow.u
            vol.water_mflow.u ~ PLR * core.mWat_flow.u
        end

        preHea.Q_flow.u ~ Qflow.u
        connect(preHea.port, vol.heatport)

        QSen_flow.u ~ core.SHR.u * Qflow.u
        QLat_flow.u ~ (1 - core.SHR.u) * Qflow.u
        P.u ~ -core.EIR.u * Qflow.u
    end
end
