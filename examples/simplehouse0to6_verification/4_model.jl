include(joinpath(@__DIR__, "3_model.jl"))

# Shared heating-loop backbone for the derived SimpleHouse4-6 cases.
# This core keeps the hydronic topology without imposing a fixed control input.
@mtkmodel SimpleHouse4Core begin
    @parameters begin
        QHea_nominal    = 3000.0
        mWat_nominal    = 0.1
        T_a_rad_nominal = 333.15
        T_b_rad_nominal = 313.15
        T_air_nominal   = 293.15
    end

    @extend weaBus, TOut, walRes, walCap, gaiWin, win, conRes, zon, RA_tap,
            AWall, dWall, kWall, rhoWall, cpWall, AWin, VZone, hWall =
            base = SimpleHouse3Core(df_weather = df_weather)

    @components begin
        heaWat = WaterHeaterCooler_Q(;
            if_vol_steady_state = true,
            Qflow_nominal = QHea_nominal,
            mflow_nominal = mWat_nominal,
        )
        rad = Radiator(;
            N_elements    = 5,
            Qflow_nominal = QHea_nominal,
            T_a_nominal   = T_a_rad_nominal,
            T_b_nominal   = T_b_rad_nominal,
            T_air_nominal = T_air_nominal,
        )
        pum    = Pump_mflow()
        bouWat = WaterBoundaryNode2_p()
    end

    @equations begin
        connect(rad.port_b,    pum.port_a)
        connect(pum.port_b,    bouWat.port_a)
        connect(bouWat.port_b, heaWat.port_a)
        connect(heaWat.port_b, rad.port_a)

        connect(rad.heatPortCon, zon.heatport)
        connect(rad.heatPortRad, walCap.port)
    end
end

# SimpleHouse4: SimpleHouse4Core + infiltration + constant-on heating.
@mtkmodel SimpleHouse4 begin
    @parameters begin
        mAir_flow_infil = 1e-5
    end

    @extend weaBus, TOut, walRes, walCap, gaiWin, win, conRes, zon, RA_tap,
            heaWat, rad, pum, bouWat, AWall, dWall, kWall, rhoWall,
            cpWall, AWin, VZone, hWall, QHea_nominal,
            mWat_nominal, T_a_rad_nominal, T_b_rad_nominal, T_air_nominal =
            base = SimpleHouse4Core(df_weather = df_weather)

    @components begin
        source = AirBoundary(role = :source, mode = :mflow, use_w_in = true)
        sink   = AirBoundary(role = :sink, mode = :pressure)
    end

    @equations begin
        source.mflow.u ~ mAir_flow_infil
        source.T_in.u  ~ zon.T_mix
        source.w_in.u  ~ zon.w_mix
        sink.p.u       ~ 101325.0
        connect(source.port, zon.airport_1)
        connect(zon.airport_2, sink.port)

        heaWat.u.u ~ 1.0
        pum.mflow_set.u ~ mWat_nominal
    end
end

