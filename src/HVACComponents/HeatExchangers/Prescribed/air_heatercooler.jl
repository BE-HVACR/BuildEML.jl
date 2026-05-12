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
