include(joinpath(@__DIR__, "2_heatingonly_model.jl"))

@mtkmodel CoolingVentLoopCore begin
    @parameters begin
        mAir_flow_nominal = 0.2
        dpAir_nominal = 200.0
        conDam_k = 1.0
        Tset_cool = unit_C2K(24.0)
        Tsup_cool = unit_C2K(20.0)
    end

    @components begin
        souSup = AirBoundary(
            role = :source,
            mode = :pressure,
            use_w_in = false,
            p_par = 101325.0,
            w_par = 0.008
        )
        sink = AirBoundary(role = :sink, mode = :pressure)
        fan = Fan_dp(
            if_dynamic_eff = false,
            mflow_nominal = mAir_flow_nominal,
            dp_nominal = dpAir_nominal
        )
        dam = DamperExponential(
            mflow_nominal = mAir_flow_nominal,
            dpDamper_nominal = 10.0,
            dpFixed_nominal = dpAir_nominal - 10.0,
            use_exponential = false
        )
        senRex = Constant_HX(
            if_latent = false,
            if_dp_const = false,
            if_dp_quadratic = false,
            eff_sensible = 0.85,
            mflow_nominal = mAir_flow_nominal,
            dp_nominal = dpAir_nominal
        )
        cooAir = AirSensibleCooler_T()
        coolingsource = ASHP_LiftCOP_Power(T_hot_offset = 8.0, T_cold_offset = -2.0)
        conDam = LimPID(k = conDam_k, u_min = 0.25, u_max = 1.0)
        set_cool = Constant(k = Tset_cool)
        set_sup = Constant(k = Tsup_cool)

        T_oa = RealInput()
        T_room = RealInput()
        totPower = RealOutput()
    end

    @equations begin
        connect(souSup.T_in, T_oa)
        sink.p.u ~ 101325.0

        connect(conDam.reference, T_room)
        connect(conDam.measurement, set_cool.output)
        connect(conDam.ctr_output, dam.y)

        connect(cooAir.TSet, set_sup.output)
        connect(coolingsource.Q_thermal, cooAir.Qflow)
        connect(coolingsource.T_hot_ref, T_oa)
        connect(coolingsource.T_cold_ref, set_sup.output)

        fan.dp_set.u ~ dpAir_nominal

        connect(souSup.port, senRex.port_a1)
        connect(senRex.port_b1, fan.port_a)
        connect(fan.port_b, cooAir.port_a)
        connect(cooAir.port_b, dam.port_a)
        connect(senRex.port_b2, sink.port)

        totPower.u ~ fan.power.u + coolingsource.power.u
    end
end

@mtkmodel SimpleHouseTotal begin
    @structural_parameters begin
        df_weather
    end

    @parameters begin
        AWall = 100.0
        VZone = AWall * 3.0
        mAir_flow_nominal = VZone * 2.0 * 1.2 / 3600.0
        dpAir_nominal = 200.0
        mflow_off_frac = 1e-4
    end

    @components begin
        house = SimpleHouseCore(
            df_weather = df_weather,
            if_infiltration = true,
            n_airports = 4,
            port_role = [:inlet, :outlet, :inlet, :outlet],
        )
        loop_hea = HeatingLoopCore()
        ctrller_hea = HeatingLoopHysteresisController()
        loop_cool = CoolingVentLoopCore(
            mAir_flow_nominal = mAir_flow_nominal,
            dpAir_nominal = dpAir_nominal,
        )
        totPower = RealOutput()
    end

    @equations begin
        connect(loop_hea.rad.heatPortCon, house.zon.heatport)
        connect(loop_hea.rad.heatPortRad, house.walCap.port)
        ctrller_hea.T_room.u ~ house.RA_tap.T.u
        loop_hea.hea.u.u ~ ctrller_hea.y
        loop_hea.pump.mflow_set.u ~
            (mflow_off_frac + (1.0 - mflow_off_frac) * ctrller_hea.y) * loop_hea.pump.mflow_nominal

        connect(loop_cool.T_oa, house.weaBus.TDryBul)
        connect(loop_cool.T_room, house.RA_tap.T)

        connect(loop_cool.dam.port_b, house.zon.airport_3)
        connect(house.zon.airport_4, loop_cool.senRex.port_a2)

        totPower.u ~ loop_hea.totPower.u + loop_cool.totPower.u
    end
end

