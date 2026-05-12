include(joinpath(@__DIR__, "0_model.jl"))

# SimpleHouse1: outdoor temperature → wall resistor → wall capacitor.
# No zone or internal air. Equivalent to Modelica's SimpleHouse1.
@mtkmodel SimpleHouse1 begin
    @parameters begin
        AWall   = 100.0
        dWall   = 0.25
        kWall   = 0.04
        rhoWall = 2000.0
        cpWall  = 1000.0
    end

    @extend SimpleHouse0(df_weather = df_weather)

    @components begin
        walRes = ThermalResistor(R = dWall / (AWall * kWall))
        walCap = HeatCapacitor(C = AWall * dWall * cpWall * rhoWall)
    end

    @equations begin
        connect(TOut.port, walRes.port_a)
        connect(walRes.port_b, walCap.port)
    end
end
