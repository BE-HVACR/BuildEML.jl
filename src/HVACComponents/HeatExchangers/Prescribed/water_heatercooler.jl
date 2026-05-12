"""
Water heater/cooler with normalized power input.

`u.u` ∈ [-1, 1] scales `Qflow_nominal`; the heat is applied to an internal mixing
volume sized by `mflow_nominal × tau`.
Ports: `port_a`/`port_b` (water), `u` (RealInput, normalized), `Qflow` (RealOutput [W]).
"""
@mtkmodel WaterHeaterCooler_Q begin
    @structural_parameters begin
        if_vol_steady_state::Bool = false
    end

    @parameters begin
        Qflow_nominal = 1000.0                           # [W] heat power at u=1
        mflow_nominal = 0.01                             # [kg/s]
        tau           = 30.0                             # [s]
        V_tank        = mflow_nominal * tau / rho_water  # [m^3]
        k_dp          = 0.0                              # [-] pressure-drop coefficient
    end

    @components begin
        port_a = WaterPort()
        port_b = WaterPort()

        vol = WaterMixingVolumeNode2(V = V_tank, if_steady_state = if_vol_steady_state)
        dp  = WaterPressureDrop(k = k_dp)

        u      = RealInput()
        Qflow  = RealOutput()

        preHea = PrescribedHeatFlow()
        gai    = Gain(k = Qflow_nominal)
    end

    @equations begin
        connect(port_a,    dp.port_a)
        connect(dp.port_b, vol.port1)
        connect(vol.port2, port_b)

        gai.u           ~ u.u
        preHea.Q_flow.u ~ gai.y
        preHea.port.T   ~ vol.T_mix
        connect(preHea.port, vol.heatport)

        Qflow.u ~ gai.y
    end
end
