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
    end
end
