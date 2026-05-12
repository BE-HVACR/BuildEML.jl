include(joinpath(@__DIR__, "4_model.jl"))

# SimpleHouse5: SimpleHouse5Core + infiltration.

@mtkmodel HeatingHysteresisCtrl5 begin
    @parameters begin
        T_low  = unit_C2K(21.0)
        T_high = unit_C2K(23.0)
    end

    @components begin
        T_room = RealInput()
        output = RealOutput()
    end

    @variables begin
        y(t)
    end

    @equations begin
        D(y) ~ ifelse(T_room.u <= T_low, 1.0, ifelse(T_room.u >= T_high, 0.0, y)) - y
        output.u ~ y
    end
end

@mtkmodel SimpleHouse5Core begin
    @parameters begin
        mflow_off_frac = 1e-4
    end

    @extend weaBus, TOut, walRes, walCap, gaiWin, win, conRes, zon, RA_tap,
            heaWat, rad, pum, bouWat, AWall, dWall, kWall, rhoWall,
            cpWall, AWin, VZone, hWall, QHea_nominal,
            mWat_nominal, T_a_rad_nominal, T_b_rad_nominal, T_air_nominal =
            base = SimpleHouse4Core(df_weather = df_weather)

    @components begin
        ctrllerHea = HeatingHysteresisCtrl5()
    end

    @equations begin
        connect(ctrllerHea.T_room, RA_tap.T)
        connect(heaWat.u, ctrllerHea.output)
        pum.mflow_set.u ~ (mflow_off_frac + (1.0 - mflow_off_frac) * ctrllerHea.output.u) * mWat_nominal
    end
end

@mtkmodel SimpleHouse5 begin
    @parameters begin
        mAir_flow_infil = 1e-5
    end

    @extend weaBus, TOut, walRes, walCap, gaiWin, win, conRes, zon, RA_tap,
            heaWat, rad, pum, bouWat, ctrllerHea, AWall, dWall, kWall, rhoWall,
            cpWall, AWin, VZone, hWall, QHea_nominal, mWat_nominal,
            T_a_rad_nominal, T_b_rad_nominal, T_air_nominal, mflow_off_frac =
            base = SimpleHouse5Core(df_weather = df_weather)

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
    end
end

