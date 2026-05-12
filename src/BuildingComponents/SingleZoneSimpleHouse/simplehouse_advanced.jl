"""
Single-zone building core with full moisture tracking, occupancy/lighting schedules, and infiltration.

Moisture sources always present: occupant latent load (via `latentToMoisture`) and infiltration (via air port humidity).
`if_extra_moisture_input=true` exposes an additional `extra_water_mflow` port for auxiliary sources (e.g. humidifier).
"""
@mtkmodel SimpleHouseCoreAdvanced begin
    @structural_parameters begin
        df_weather
        df_schedule
        if_infiltration::Bool = true
        if_extra_moisture_input::Bool = false
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
        h_fg_occ = 2.5e6
        ACH = 0.2
        mAir_flow_infil = 1.2 * VZone * ACH / 3600
    end

    @components begin
        zon = AirMixingVolumeNodeN(
            if_moisture_input = true,
            N = n_airports,
            port_role = port_role,
            if_steady_state = false,
            V = VZone,
            Qflow_const = 0.0,
        )
        weaBus = WeatherBus(df_weather, interp_method = DataInterpolations.AkimaInterpolation)
        sch = ScheduleBus(df_schedule)
        conRes = ThermalResistor(R = R_con)
        walRes = ThermalResistor(R = R_wal)
        walCap = HeatCapacitor(C = C_wal)
        TOut = PrescribedTemperature()
        gaiWin = Gain(k = AWin * gWin)
        win = PrescribedHeatFlow()
        occSen = PrescribedHeatFlow()
        occLat = PrescribedHeatFlow()
        lightCon = PrescribedHeatFlow()
        lightRad = PrescribedHeatFlow()
        latentToMoisture = Gain(k = 1.0 / h_fg_occ)
        senT = TemperatureSensor()
        RA_tap = AirStateTapCore()
        source = AirBoundary(role = :source, mode = :mflow, use_w_in = true)
        if if_extra_moisture_input
            extra_water_mflow = RealInput()
        end
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

        connect(sch.QOccSen, occSen.Q_flow)
        connect(sch.QOccLat, occLat.Q_flow)
        connect(sch.QLightCon, lightCon.Q_flow)
        connect(sch.QLightRad, lightRad.Q_flow)
        connect(sch.QOccLat, latentToMoisture.input)
        if if_extra_moisture_input
            zon.water_mflow.u ~ latentToMoisture.output.u + extra_water_mflow.u
        else
            zon.water_mflow.u ~ latentToMoisture.output.u
        end

        connect(TOut.port, walRes.port_a)
        connect(walRes.port_b, walCap.port)
        connect(zon.heatport, conRes.port_a)
        connect(conRes.port_b, walCap.port)
        connect(win.port, walCap.port)
        connect(lightRad.port, walCap.port)
        connect(occSen.port, zon.heatport)
        connect(occLat.port, zon.heatport)
        connect(lightCon.port, zon.heatport)
        connect(senT.port, zon.heatport)

        RA_tap.T_probe.u ~ zon.T_mix
        RA_tap.w_probe.u ~ zon.w_mix

        if if_infiltration
            source.mflow.u ~ mAir_flow_infil
            connect(source.T_in, weaBus.TDryBul)
            connect(source.w_in, weaBus.HumRat)
            if n_airports == 2
                sink.p.u ~ 101325.0
            else
                sink.mflow.u ~ mAir_flow_infil
            end
            connect(source.port, zon.airport_1)
            connect(zon.airport_2, sink.port)
        end
    end
end

