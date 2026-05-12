include(joinpath(@__DIR__, "2_model.jl"))

# Shared zone-air backbone for the derived SimpleHouse3-6 cases.
# This core adds the zone, wall convection, and state taps without prescribing an air-side topology.
@mtkmodel SimpleHouse3Core begin
    @parameters begin
        VZone = 8.0 * 8.0 * 3.0
        hWall = 2.0
    end

    @extend weaBus, TOut, walRes, walCap, gaiWin, win, AWall, dWall, kWall, rhoWall, cpWall, AWin =
            base = SimpleHouse2(df_weather = df_weather)

    @components begin
        conRes = ThermalResistor(R = 1.0 / (hWall * AWall))
        zon    = AirMixingVolumeNodeN(N = 2, if_steady_state = false, V = VZone, Qflow_const = 0.0)
        RA_tap = AirStateTapCore()
    end

    @equations begin
        connect(walCap.port, conRes.port_b)
        connect(conRes.port_a, zon.heatport)

        RA_tap.T_probe.u ~ zon.T_mix
        RA_tap.w_probe.u ~ zon.w_mix
    end
end

# SimpleHouse3: SimpleHouse3Core + tiny infiltration to keep the air port numerically stable.
@mtkmodel SimpleHouse3 begin
    @parameters begin
        mAir_flow_infil = 1e-5
    end

    @extend weaBus, TOut, walRes, walCap, gaiWin, win, conRes, zon, RA_tap,
            AWall, dWall, kWall, rhoWall, cpWall, AWin, VZone, hWall =
            base = SimpleHouse3Core(df_weather = df_weather)

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

