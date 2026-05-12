@mtkmodel HeatingLoopCore begin
    @parameters begin
        heating_capacity = 700.0
        Tdrop_radiator_nominal = 10.0
        mflow_set_par = heating_capacity / Tdrop_radiator_nominal / cp_water
    end

    @components begin
        bou = WaterBoundaryNode2_p(use_p_in = false)

        pump = Pump_mflow(
            mflow_nominal = mflow_set_par,
            dp_nominal    = 5000.0,
        )

        hea  = WaterHeaterCooler_Q(
            if_vol_steady_state = true,
            Qflow_nominal = heating_capacity,
            mflow_nominal = mflow_set_par,
        )

        rad  = Radiator(
            N_elements    = 5,
            Qflow_nominal = heating_capacity,
            T_a_nominal   = unit_C2K(50.0),
            T_b_nominal   = unit_C2K(50.0) - Tdrop_radiator_nominal,
            T_air_nominal = 293.15,
            dp_nominal    = 5000.0
        )

        heatingsrc_power = ConstantCOP_Power(COP_const = 1.0)
        totPower = RealOutput()
    end

    @equations begin
        connect(bou.port_b, hea.port_a)
        connect(hea.port_b,   rad.port_a)
        connect(rad.port_b,   pump.port_a)
        connect(pump.port_b,  bou.port_a)

        connect(heatingsrc_power.Q_thermal, hea.Qflow)
        totPower.u ~ heatingsrc_power.power.u + pump.power.u
    end
end;

@mtkmodel HeatingLoopHysteresisController begin
    @parameters begin
        T_low  = unit_C2K(20.0)   # [K]
        T_high = unit_C2K(22.0)   # [K]
    end

    @components begin
        T_room = RealInput()
        output  = RealOutput()
    end

    @variables begin
        y(t)
    end

    @equations begin
        D(y) ~ ifelse(T_room.u <= T_low, 1.0, ifelse(T_room.u >= T_high, 0.0, y)) - y
        output.u ~ y
    end
end

@mtkmodel SimpleHouseHeatingOnly begin
    @structural_parameters begin
        df_weather
    end

    @parameters begin
        Tset_room = unit_C2K(22.0)
        ΔT_band   = 1.0
        mflow_off_frac = 1e-4
    end

    @components begin
        house = SimpleHouseCore(
            df_weather = df_weather,
            n_airports = 2,
            port_role = [:inlet, :outlet],
        )
        loop  = HeatingLoopCore()
        ctrller_hea = HeatingLoopHysteresisController()
    end

    @equations begin
        connect(loop.rad.heatPortCon, house.zon.heatport)
        connect(loop.rad.heatPortRad, house.walCap.port)

        ctrller_hea.T_room.u ~ house.RA_tap.T.u

        loop.hea.u.u           ~ ctrller_hea.y
        loop.pump.mflow_set.u  ~ (mflow_off_frac + (1.0 - mflow_off_frac) * ctrller_hea.y) * loop.pump.mflow_nominal
    end
end
