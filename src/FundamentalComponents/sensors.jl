@mtkmodel AirStateTapCore begin
    @components begin
        T_probe  = RealInput()
        w_probe  = RealInput()
        T = RealOutput()
        w = RealOutput()

    end
    @equations begin
        connect(T_probe, T)
        connect(w_probe, w)
    end
end
