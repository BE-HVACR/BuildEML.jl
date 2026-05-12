"""
Single-zone building core with wall conduction, solar gain, and optional infiltration.
Zone humidity (w_mix) is always tracked. `use_w_in` controls whether infiltrating air
carries outdoor humidity (true) or a fixed default w_par (false).
"""
@mtkmodel SimpleHouseCore begin
    @structural_parameters begin
        df_weather
        if_infiltration::Bool = true
        use_w_in::Bool = true
        n_airports::Int = 2
        port_role = [:inlet, :outlet]
    end

    @parameters begin
        AWall = 100.0
        AWin = 5.0
        gWin = 0.3
        VZone = AWall * 3.0
        d_wal = 0.25
        k_wal = 0.04
        rho_wal = 2000.0
        cp_wal = 1000.0
        h_in_wal = 2.0
        R_con = 1.0 / (h_in_wal * AWall)
        R_wal = d_wal / (AWall * k_wal)
        C_wal = AWall * d_wal * rho_wal * cp_wal
        ACH = 1e-3
        mAir_flow_infil = 1.2 * VZone * ACH / 3600
    end

    @components begin
        zon = AirMixingVolumeNodeN(
            N = n_airports,
            port_role = port_role,
            if_steady_state = false,
            V = VZone,
            Qflow_const = 0.0,
        )
        weaBus = WeatherBus(df_weather, interp_method = DataInterpolations.AkimaInterpolation)
        conRes = ThermalResistor(R = R_con)
        walRes = ThermalResistor(R = R_wal)
        walCap = HeatCapacitor(C = C_wal)
        TOut = PrescribedTemperature()
        gaiWin = Gain(k = AWin * gWin)
        win = PrescribedHeatFlow()
        senT = TemperatureSensor()
        RA_tap = AirStateTapCore()
        source = AirBoundary(role = :source, mode = :mflow, use_w_in = use_w_in)
        if n_airports == 2
            sink = AirBoundary(role = :sink, mode = :pressure)
        else
            sink = AirBoundary(role = :sink, mode = :mflow)
        end
    end

    @equations begin
        connect(weaBus.TDryBul, TOut.T)
        connect(weaBus.HGloHor, gaiWin.input)
        connect(gaiWin.output, win.Q_flow)

        connect(TOut.port, walRes.port_a)
        connect(walRes.port_b, walCap.port)
        connect(zon.heatport, conRes.port_a)
        connect(conRes.port_b, walCap.port)
        connect(win.port, walCap.port)
        connect(senT.port, zon.heatport)

        RA_tap.T_probe.u ~ zon.T_mix
        RA_tap.w_probe.u ~ zon.w_mix

        if if_infiltration
            source.mflow.u ~ mAir_flow_infil
            connect(source.T_in, weaBus.TDryBul)
            if n_airports == 2
                sink.p.u ~ 101325.0
            else
                sink.mflow.u ~ mAir_flow_infil
            end
            connect(source.port, zon.airport_1)
            connect(zon.airport_2, sink.port)
        end
        if if_infiltration && use_w_in
            connect(source.w_in, weaBus.HumRat)
        end
    end
end

