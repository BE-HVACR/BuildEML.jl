"""Quadratic water pressure-drop resistance: dp = k · ṁ · |ṁ|."""

@mtkmodel WaterPressureDrop begin
    @parameters begin
        k = 1.0
    end
    @components begin
        port_a = WaterPort()
        port_b = WaterPort()
    end
    @variables begin
        mflow(t)
        dp(t)
    end
    @equations begin
        mflow ~ port_a.mflow
        port_a.mflow + port_b.mflow ~ 0

        port_b.p ~ port_a.p - dp
        dp ~ k * mflow * smooth_abs(mflow)

        port_b.T_ofo ~ instream(port_a.T_ofo)
        port_a.T_ofo ~ port_a.T_guess
    end
end


"""
Air flow resistance with a square-root law and a polynomial transition at low flow.

The total (moist-air) mass flow `mflow * (1 + w_in)` follows MBL's
`basicFlowFunction_dp`, as in `Buildings.Fluid.FixedResistances.PressureDrop`
with `linearized = false`: `mflow_nominal` at `dp_nominal`, and the transition
below `deltaM * mflow_nominal` (same default 0.3). `mflow_nominal` is a
moist-air mass flow rate.
"""
@mtkmodel AirPressureDrop begin
    @parameters begin
        mflow_nominal = 1.0     # [kg/s] moist air
        dp_nominal = 100.0      # [Pa]
        deltaM = 0.3            # [-]
        k = mflow_nominal / sqrt(dp_nominal)
        m_turb = deltaM * mflow_nominal
        dp_turb = (m_turb / k)^2
    end
    @components begin
        port_a = AirPort()
        port_b = AirPort()
    end
    @variables begin
        mflow(t)                    # [kg_da/s]
        dp(t)                       # [Pa]
        T_in(t), [guess = 293.15]   # [K]
        w_in(t), [guess = 0.008]    # [kg_w/kg_da]
    end
    @equations begin
        mflow ~ port_a.mflow
        port_a.mflow + port_b.mflow ~ 0
        dp ~ port_a.p - port_b.p

        T_in ~ instream(port_a.T_ofo)
        w_in ~ instream(port_a.w_ofo)
        mflow * (1 + w_in) ~ ifelse(abs(dp) > dp_turb,
            sign(dp) * k * sqrt(abs(dp)),
            (1.40625 + (0.15625 * (dp / dp_turb)^2 - 0.5625) * (dp / dp_turb)^2) *
                m_turb * (dp / dp_turb))

        port_b.T_ofo ~ T_in
        port_b.w_ofo ~ w_in
        port_a.T_ofo ~ port_a.T_guess
        port_a.w_ofo ~ port_a.w_guess
    end
end
