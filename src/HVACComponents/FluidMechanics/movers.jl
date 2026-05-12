"""
Ideal constant-pressure fan (FlowControlled_dp-type).

Ports: `port_a` (inlet), `port_b` (outlet), `dp_set` (pressure rise [Pa]), `power` (electrical power [W]).
Air state (T, w) passes through unchanged; pressure rise follows `dp_set.u`.

Power is computed from dp × Vflow and efficiency. If `if_dynamic_eff = true`,
`variablefan_type` selects the PLR-based curve: "VSD" (default), "InletVane", or "OutletDamper".
"""
@mtkmodel Fan_dp begin
    @structural_parameters begin
        if_dynamic_eff::Bool = false
        variablefan_type::String = "VSD"
    end

    @parameters begin
        rho_nominal           = 1.2            # [kg/m^3]
        dp_nominal            = 200.0          # [Pa]
        mflow_nominal         = 0.2            # [kg/s]
        eff_hydraulic_nominal = 0.7            # [-]
        eff_motor_nominal     = 0.7            # [-]
        eff_tot_nominal       = 0.7 * 0.7     # [-]
        power_nominal = 200.0 * 0.2 / (1.2 * (0.7 * 0.7))  # [W]
    end

    @components begin
        port_a = AirPort()
        port_b = AirPort()
        dp_set = RealInput()    # [Pa]
        power  = RealOutput()   # [W]
    end

    @variables begin
        mflow(t)   # [kg/s]
        Vflow(t)   # [m^3/s]
        dp_fan(t)  # [Pa]
        power_(t)  # [W]
        PLR(t)
    end

    @equations begin
        port_b.T_ofo ~ instream(port_a.T_ofo)
        port_b.w_ofo ~ instream(port_a.w_ofo)

        dp_fan   ~ dp_set.u
        port_b.p ~ dp_fan + port_a.p

        port_a.mflow + port_b.mflow ~ 0
        mflow ~ port_a.mflow

        Vflow ~ mflow / rho_nominal

        if !(if_dynamic_eff)
            power_ ~ dp_fan * Vflow / eff_tot_nominal
        else
            PLR ~ smooth_clamp01(mflow / mflow_nominal)
            if variablefan_type == "VSD"
                power_ ~ power_nominal * (0.00153 + 0.0052 * PLR + 1.1086 * PLR^2 - 0.1164 * PLR^3)
            elseif variablefan_type == "InletVane"
                power_ ~ power_nominal * (0.351 + 0.308 * PLR - 0.541 * PLR^2 + 0.872 * PLR^3)
            elseif variablefan_type == "OutletDamper"
                power_ ~ power_nominal * (0.371 + 0.973 * PLR - 0.3342 * PLR^2)
            end
        end

        power.u ~ power_
    end
end


"""
Ideal mass-flow-controlled water pump.

`mflow_set.u` prescribes the mass flow; pressure rise adapts to the connected network.
Positive flow direction is port_a → port_b.
Power is computed from hydraulic work and efficiency (constant or PLR-based via `if_dynamic_eff`).
"""
@mtkmodel Pump_mflow begin
    @structural_parameters begin
        if_dynamic_eff::Bool = false
    end

    @parameters begin
        mflow_nominal         = 0.1    # [kg/s]
        dp_nominal            = 5000.0 # [Pa]
        eff_hydraulic_nominal = 0.7    # [-]
        eff_motor_nominal     = 0.7    # [-]
        eff_tot_nominal       = eff_hydraulic_nominal * eff_motor_nominal  # [-]
        power_nominal         = mflow_nominal * dp_nominal / (rho_water * eff_tot_nominal)  # [W]
    end

    @components begin
        port_a    = WaterPort()
        port_b    = WaterPort()
        mflow_set = RealInput()    # [kg/s]
        power     = RealOutput()   # [W]
    end

    @variables begin
        mflow(t)          # [kg/s]
        dp_pump(t)        # [Pa]  port_b.p - port_a.p
        Vflow(t)          # [m^3/s]
        PLR(t)            # [-]   part-load ratio
        eff_hydraulic(t)  # [-]
        eff_motor(t)      # [-]
        power_(t)         # [W]
    end

    @equations begin
        port_b.T_ofo ~ instream(port_a.T_ofo)

        mflow ~ mflow_set.u
        port_a.mflow ~  mflow
        port_b.mflow ~ -mflow

        dp_pump ~ port_b.p - port_a.p
        Vflow   ~ mflow / rho_water

        if !(if_dynamic_eff)
            power_ ~ dp_pump * Vflow / eff_tot_nominal
        else
            PLR ~ smooth_clamp01(mflow / mflow_nominal)
            eff_hydraulic ~ eff_hydraulic_nominal * (0.25 + 0.75 * PLR - 0.10 * PLR^2)
            eff_motor     ~ eff_motor_nominal     * (0.35 + 0.65 * PLR - 0.05 * PLR^2)
            power_ ~ dp_pump * Vflow / (eff_hydraulic * eff_motor)
        end

        power.u ~ power_
    end
end
