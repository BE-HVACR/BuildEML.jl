"""
Sensible cooler with prescribed outlet temperature.

- One air stream: port_a -> port_b
- Only sensible heat exchange: humidity ratio passes through unchanged.
- Cooling-only: Qflow ≤ 0. The requested cooling from TSet is limited to [Qflow_lb, 0].
"""
@mtkmodel AirSensibleCooler_T begin
    @parameters begin
        Qflow_lb  = -1e9    # [W]
        mcp_min   = 1e-4    # [W/K] regularization lower bound
    end

    @components begin
        port_a = AirPort()
        port_b = AirPort()

        TSet  = RealInput()    # [K]
        Qflow = RealOutput()   # [W] cooling < 0
    end

    @variables begin
        mflow(t)                        # [kg/s]
        T_in(t),  [guess = 293.15]
        T_out(t), [guess = 293.15]
        mcp_eff(t)                      # [W/K]
        Qflow_req(t)                    # [W]
        Qflow_actual(t)                 # [W]
    end

    @equations begin
        port_a.mflow + port_b.mflow ~ 0
        mflow ~ port_a.mflow

        port_b.p ~ port_a.p

        T_in ~ instream(port_a.T_ofo)

        Qflow_req ~ mflow * cp_da * (TSet.u - T_in)

        # Cooling-only with capacity limit:
        #   Qflow_actual ∈ [Qflow_lb, 0]
        #   (Q>0 would heat; we force max 0)
        Qflow_actual ~ smooth_clamp(Qflow_req, Qflow_lb, 0.0)

        # regularize near zero flow to avoid singular initialization
        mcp_eff ~ smooth_max(smooth_abs(mflow * cp_da), mcp_min)
        T_out   ~ T_in + Qflow_actual / mcp_eff

        port_b.T_ofo ~ T_out
        port_b.w_ofo ~ instream(port_a.w_ofo)

        Qflow.u ~ Qflow_actual
    end
end

"""
Sensible heater with prescribed outlet temperature.

- One air stream: `port_a -> port_b`
- Only sensible heat exchange: humidity ratio passes through unchanged.
- Heating-only: `Qflow >= 0`. The requested heating from `TSet` is limited to `[0, Qflow_ub]`.
"""
@mtkmodel AirSensibleHeater_T begin
    @parameters begin
        Qflow_ub = 1e9
        mcp_min  = 1e-4
    end

    @components begin
        port_a = AirPort()
        port_b = AirPort()

        TSet  = RealInput()
        Qflow = RealOutput()
    end

    @variables begin
        mflow(t)
        T_in(t),  [guess = 293.15]
        T_out(t), [guess = 293.15]
        mcp_eff(t)
        Qflow_req(t)
        Qflow_actual(t)
    end

    @equations begin
        port_a.mflow + port_b.mflow ~ 0
        mflow ~ port_a.mflow

        port_b.p ~ port_a.p
        T_in ~ instream(port_a.T_ofo)

        Qflow_req    ~ mflow * cp_da * (TSet.u - T_in)
        Qflow_actual ~ smooth_clamp(Qflow_req, 0.0, Qflow_ub)

        mcp_eff ~ smooth_max(smooth_abs(mflow * cp_da), mcp_min)
        T_out   ~ T_in + Qflow_actual / mcp_eff

        port_b.T_ofo ~ T_out
        port_b.w_ofo ~ instream(port_a.w_ofo)

        Qflow.u ~ Qflow_actual
    end
end

"""
Wet cooling coil with prescribed apparatus dew point.

- One air stream: port_a -> port_b
- Simultaneous sensible and latent cooling.
- Outlet: T_out = min(T_coil, T_in), w_out = min(w_sat(T_out), w_in) — dehumidification only, no humidification.
- T_coil = max(T_coil_sp.u, T_dp_min): coil cannot cool below its physical dew point limit.
- Qflow is total enthalpy removal (sensible + latent), always ≤ 0.
"""
@mtkmodel AirWetCooler begin
    @parameters begin
        T_dp_min = 277.15   # [K] minimum achievable apparatus dew point (e.g. 4°C)
        mcp_min  = 1e-4     # [W/K] regularization lower bound
    end

    @components begin
        port_a    = AirPort()
        port_b    = AirPort()
        T_coil_sp = RealInput()    # [K] apparatus dew point setpoint
        Qflow     = RealOutput()   # [W] total heat removal (≤ 0)
    end

    @variables begin
        mflow(t)
        T_in(t),   [guess = 293.15]
        w_in(t),   [guess = 0.010]
        T_coil(t), [guess = 285.15]
        T_out(t),  [guess = 285.15]
        w_out(t),  [guess = 0.009]
        Qflow_total(t)
    end

    @equations begin
        port_a.mflow + port_b.mflow ~ 0
        mflow ~ port_a.mflow

        port_b.p ~ port_a.p

        T_in ~ instream(port_a.T_ofo)
        w_in ~ instream(port_a.w_ofo)

        T_coil ~ smooth_max(T_coil_sp.u, T_dp_min)

        # Cooling-only: coil cannot be warmer than the inlet air
        T_out ~ smooth_min(T_coil, T_in)
        w_out ~ smooth_min(w_sat(T_out), w_in)

        Qflow_total ~ mflow * (h_Tw(T_out, w_out) - h_Tw(T_in, w_in))

        port_b.T_ofo ~ T_out
        port_b.w_ofo ~ w_out

        Qflow.u ~ Qflow_total
    end
end

"""
Control block: selects the colder of two competing coil temperature demands.

- T_coil_sp = smooth_min(T_hum_sp, T_sa_sp)
  - T_sa_sp colder: temperature drives the coil, no reheat needed.
  - T_hum_sp colder: humidity drives the coil, reheat back to T_sa_sp.
- T_hum_sp is output directly by the humidity PI in temperature units.
- The physical lower limit T_dp_min is enforced inside AirWetCooler, not here.
"""
@mtkmodel CoilTempSelector begin
    @components begin
        T_hum_sp  = RealInput()    # [K] humidity-driven coil temperature setpoint (from humidity PI)
        T_sa_sp   = RealInput()    # [K] supply air temperature setpoint
        T_coil_sp = RealOutput()   # [K] apparatus dew point setpoint → AirWetCooler.T_coil_sp
    end

    @equations begin
        T_coil_sp.u ~ smooth_min(T_hum_sp.u, T_sa_sp.u)
    end
end
