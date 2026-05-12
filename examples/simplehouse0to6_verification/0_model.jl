using BuildEML

epw_path = joinpath(@__DIR__, "..", "..", "src", "Disturbances", "Weather", "weatherfile",
    "USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.epw")
df_weather = ReadEPW(epw_path)

@mtkmodel SimpleHouse0 begin
    @structural_parameters begin
        df_weather
    end

    @components begin
        weaBus = WeatherBus(df_weather, interp_method = DataInterpolations.AkimaInterpolation)
        TOut   = PrescribedTemperature()
    end

    @equations begin
        connect(weaBus.TDryBul, TOut.T)
    end
end

