include(joinpath(@__DIR__, "1_model.jl"))

# SimpleHouse2: SimpleHouse1 + window solar gain (HDirNor * AWin) injected onto walCap.
@mtkmodel SimpleHouse2 begin
    @parameters begin
        AWin = 2.0
    end

    @extend weaBus, TOut, walRes, walCap, AWall, dWall, kWall, rhoWall, cpWall = base = SimpleHouse1(df_weather = df_weather)

    @components begin
        gaiWin = Gain(k = AWin)
        win    = PrescribedHeatFlow()
    end

    @equations begin
        connect(weaBus.HDirNor, gaiWin.input)
        connect(gaiWin.output, win.Q_flow)
        connect(win.port, walCap.port)
    end
end
