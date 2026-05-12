"""
Single water-to-air radiator element: extends `WaterMixingVolumeNode2` with a power-law UA
that splits heat output into convective (`heatPortCon`) and radiant (`heatPortRad`) fractions.
"""
@mtkmodel RadiatorSingleElement begin
    @parameters begin
        UA_element = 10.0      # [W/K^n]
        V_element  = 0.002     # [m^3]

        M_metal_element = 1.0     # [kg]
        cp_metal_       = 500.0   # [J/kg·K]

        n_exp_  = 1.24   # [-] heat transfer exponent
        fraRad_ = 0.35   # [-] radiant fraction
    end

    @components begin
        heatPortCon = HeatPort()
        heatPortRad = HeatPort()

        heatsource = PrescribedHeatFlow()
    end

    @extend port1, port2, T_mix, heatport = vol = WaterMixingVolumeNode2(V = V_element)

    @variables begin
        Q_con(t)
        Q_rad(t)

        Q_metal(t)

        dT_con(t)
        dT_rad(t)
    end

    @equations begin
        dT_con ~ T_mix - heatPortCon.T
        dT_rad ~ T_mix - heatPortRad.T

        Q_con ~ (1 - fraRad_) * UA_element * dT_con * smooth_abs(dT_con)^(n_exp_ - 1)
        Q_rad ~ fraRad_      * UA_element * dT_rad * smooth_abs(dT_rad)^(n_exp_ - 1)
        Q_metal ~ (M_metal_element * cp_metal_) * D(T_mix)

        heatsource.Q_flow.u ~ -(Q_con + Q_rad + Q_metal)
        connect(heatsource.port, heatport)

        heatPortCon.Q_flow ~ -Q_con
        heatPortRad.Q_flow ~ -Q_rad
    end
end


"""
N-element panel radiator with pressure drop.

Ports: `port_a`/`port_b` (water inlet/outlet), `heatPortCon` (convection), `heatPortRad` (radiation).
UA, water volume, and metal mass are derived from nominal operating conditions; `N_elements`
controls discretization granularity.
"""
@mtkmodel Radiator begin
    @structural_parameters begin
        N_elements::Int = 1
    end

    @parameters begin
        Qflow_nominal  = 700.0            # [W]
        T_a_nominal    = unit_C2K(50.0)   # [K]
        T_b_nominal    = unit_C2K(40.0)   # [K]
        T_air_nominal  = 293.15           # [K]

        dp_nominal     = 5000.0           # [Pa]

        V_water_total  = 5.8e-6 * abs(Qflow_nominal)    # [m^3]

        cp_metal      = 500.0                            # [J/kg·K]
        M_metal_total = 0.0263 * abs(Qflow_nominal)     # [kg]

        n_exp  = 1.24   # [-]
        fraRad = 0.35   # [-]

        mflow_nominal = Qflow_nominal / (cp_water * (T_a_nominal - T_b_nominal))
        k_dp = dp_nominal / (mflow_nominal^2)

        T_mean    = (T_a_nominal + T_b_nominal) / 2
        T_eff     = T_air_nominal
        dT_meanL  = T_mean - T_eff
        UA_total  = Qflow_nominal / (dT_meanL^n_exp)
        UA_ele_val = UA_total / N_elements

        V_ele_val  = V_water_total / N_elements
        M_ele_val  = M_metal_total / N_elements
    end

    @components begin
        port_a = WaterPort()
        port_b = WaterPort()

        heatPortCon = HeatPort()
        heatPortRad = HeatPort()

        elems = [
            RadiatorSingleElement(
                UA_element      = UA_ele_val,
                V_element       = V_ele_val,
                M_metal_element = M_ele_val,
                cp_metal_       = cp_metal,
                n_exp_          = n_exp,
                fraRad_         = fraRad
            ) for i in 1:N_elements
        ]

        res = WaterPressureDrop(k = k_dp)
    end

    @equations begin
        connect(port_a, res.port_a)
        connect(res.port_b, elems[1].port1)

        [connect(elems[i].port2, elems[i+1].port1) for i in 1:(N_elements-1)]...

        connect(elems[N_elements].port2, port_b)

        [connect(elems[i].heatPortCon, heatPortCon) for i in 1:N_elements]...
        [connect(elems[i].heatPortRad, heatPortRad) for i in 1:N_elements]...
    end
end
