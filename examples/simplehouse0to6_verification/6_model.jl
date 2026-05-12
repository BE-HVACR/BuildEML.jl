include(joinpath(@__DIR__, "5_model.jl"))

# SimpleHouse6: SimpleHouse5Core + free-cooling ventilation (HRV + fan + damper).
# Changes from SimpleHouse5: the effective window area increases from 2 m^2 to 6 m^2, and the infiltration branch is replaced by ventilation.
# Supply  : souVent(p=101325, T=TOut) -> hexRec(a1->b1) -> fan -> dam -> zon.port1
# Exhaust : zon.port2 -> hexRec(a2->b2) -> sinkVent(p=101325)

@mtkmodel SimpleHouse6 begin
    @parameters begin
        AWin_extra        = 4.0
        mAir_flow_nominal = 0.1
        dpAir_nominal     = 200.0
        T_low_vent        = 273.15 + 23.0
        T_high_vent       = 273.15 + 25.0
        yVent_off         = 1e-4
    end

    @variables begin
        yVent(t) = 1e-4, [irreducible = true]
    end

    @extend weaBus, TOut, walRes, walCap, gaiWin, win, conRes, zon, RA_tap,
            heaWat, rad, pum, bouWat, ctrllerHea, AWall, dWall, kWall, rhoWall,
            cpWall, VZone, hWall, QHea_nominal, mWat_nominal, T_a_rad_nominal,
            T_b_rad_nominal, T_air_nominal, mflow_off_frac =
            base = SimpleHouse5Core(df_weather = df_weather)

    @components begin
        gaiWinExtra = Gain(k = AWin_extra)
        winExtra    = PrescribedHeatFlow()
        hexRec = Constant_HX(;
            eff_sensible  = 0.85,
            if_latent     = false,
            if_dp_const   = true,
            mflow_nominal = mAir_flow_nominal,
            dp_nominal = 10.0,
        )
        fan = Fan_dp(;
            dp_nominal    = dpAir_nominal,
            mflow_nominal = mAir_flow_nominal,
        )
        dam = DamperExponential(;
            mflow_nominal    = mAir_flow_nominal,
            dpDamper_nominal = dpAir_nominal,
            dpFixed_nominal  = 0.0,
            use_exponential  = false,
        )
        souVent  = AirBoundary(role = :source, mode = :pressure, use_w_in = false, p_par = 101325.0, w_par = 0.008)
        sinkVent = AirBoundary(role = :sink, mode = :pressure)
    end

    @equations begin
        D(yVent) ~ 0.0

        connect(weaBus.HDirNor, gaiWinExtra.input)
        connect(gaiWinExtra.output, winExtra.Q_flow)
        connect(winExtra.port, walCap.port)

        connect(souVent.T_in, weaBus.TDryBul)
        sinkVent.p.u ~ 101325.0

        connect(souVent.port,   hexRec.port_a1)
        connect(hexRec.port_b1, fan.port_a)
        connect(fan.port_b,     dam.port_a)
        connect(dam.port_b,     zon.airport_1)

        connect(zon.airport_2,      hexRec.port_a2)
        connect(hexRec.port_b2, sinkVent.port)

        fan.dp_set.u ~ yVent * dpAir_nominal
        dam.y.u ~ yVent
    end
end

