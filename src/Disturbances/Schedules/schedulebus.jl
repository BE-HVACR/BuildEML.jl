"""
    ScheduleBus(df; name=:ScheduleBus, time_col=:time, interp_method=ConstantInterpolation)

Expose a schedule table as a bus of piecewise-constant signals.

This mirrors the `WeatherBus(df)` workflow:
- first prepare a DataFrame
- then feed it into a bus block consumed by the ODE/DAE model
"""
function ScheduleBus(df::DataFrame;
                     name::Symbol = :ScheduleBus,
                     time_col::Symbol = :time,
                     interp_method = ConstantInterpolation)
    _validate_schedule_dataframe(df; time_col = time_col)

    time_data = Float64.(df[!, time_col])
    occ_ratio_data = Float64.(df[!, :occRatio])
    occ_active_data = Float64.(df[!, :occActive])
    people_count_data = Float64.(df[!, :peopleCount])
    vent_min_data = Float64.(df[!, :ventMin])
    QOccSen_data = Float64.(df[!, :QOccSen])
    QOccLat_data = Float64.(df[!, :QOccLat])
    PLight_data = Float64.(df[!, :PLight])
    QLightRad_data = Float64.(df[!, :QLightRad])
    QLightCon_data = Float64.(df[!, :QLightCon])
    THeaSet_data = Float64.(df[!, :THeaSet])
    TCooSet_data = Float64.(df[!, :TCooSet])
    RHMinSet_data = Float64.(df[!, :RHMinSet])
    RHMaxSet_data = Float64.(df[!, :RHMaxSet])

    # Match the weather-bus behavior by extending the right endpoint one bin.
    # This avoids end-of-year extrapolation errors when the solver evaluates
    # exactly at the final simulation time.
    if length(time_data) >= 2
        dt = time_data[2] - time_data[1]
        time_data = vcat(time_data, time_data[end] + dt)
        occ_ratio_data = vcat(occ_ratio_data, occ_ratio_data[end])
        occ_active_data = vcat(occ_active_data, occ_active_data[end])
        people_count_data = vcat(people_count_data, people_count_data[end])
        vent_min_data = vcat(vent_min_data, vent_min_data[end])
        QOccSen_data = vcat(QOccSen_data, QOccSen_data[end])
        QOccLat_data = vcat(QOccLat_data, QOccLat_data[end])
        PLight_data = vcat(PLight_data, PLight_data[end])
        QLightRad_data = vcat(QLightRad_data, QLightRad_data[end])
        QLightCon_data = vcat(QLightCon_data, QLightCon_data[end])
        THeaSet_data = vcat(THeaSet_data, THeaSet_data[end])
        TCooSet_data = vcat(TCooSet_data, TCooSet_data[end])
        RHMinSet_data = vcat(RHMinSet_data, RHMinSet_data[end])
        RHMaxSet_data = vcat(RHMaxSet_data, RHMaxSet_data[end])
    end

    @named clk = ContinuousClock()

    @named itp_occ_ratio = Interpolation(interp_method, occ_ratio_data, time_data)
    @named itp_occ_active = Interpolation(interp_method, occ_active_data, time_data)
    @named itp_people_count = Interpolation(interp_method, people_count_data, time_data)
    @named itp_vent_min = Interpolation(interp_method, vent_min_data, time_data)
    @named itp_QOccSen = Interpolation(interp_method, QOccSen_data, time_data)
    @named itp_QOccLat = Interpolation(interp_method, QOccLat_data, time_data)
    @named itp_PLight = Interpolation(interp_method, PLight_data, time_data)
    @named itp_QLightRad = Interpolation(interp_method, QLightRad_data, time_data)
    @named itp_QLightCon = Interpolation(interp_method, QLightCon_data, time_data)
    @named itp_THeaSet = Interpolation(interp_method, THeaSet_data, time_data)
    @named itp_TCooSet = Interpolation(interp_method, TCooSet_data, time_data)
    @named itp_RHMinSet = Interpolation(interp_method, RHMinSet_data, time_data)
    @named itp_RHMaxSet = Interpolation(interp_method, RHMaxSet_data, time_data)

    @named occ_ratio = RealOutput()
    @named occ_active = RealOutput()
    @named people_count = RealOutput()
    @named vent_min = RealOutput()
    @named QOccSen = RealOutput()
    @named QOccLat = RealOutput()
    @named PLight = RealOutput()
    @named QLightRad = RealOutput()
    @named QLightCon = RealOutput()
    @named THeaSet = RealOutput()
    @named TCooSet = RealOutput()
    @named RHMinSet = RealOutput()
    @named RHMaxSet = RealOutput()

    eqs = [
        connect(clk.output, itp_occ_ratio.input)
        connect(clk.output, itp_occ_active.input)
        connect(clk.output, itp_people_count.input)
        connect(clk.output, itp_vent_min.input)
        connect(clk.output, itp_QOccSen.input)
        connect(clk.output, itp_QOccLat.input)
        connect(clk.output, itp_PLight.input)
        connect(clk.output, itp_QLightRad.input)
        connect(clk.output, itp_QLightCon.input)
        connect(clk.output, itp_THeaSet.input)
        connect(clk.output, itp_TCooSet.input)
        connect(clk.output, itp_RHMinSet.input)
        connect(clk.output, itp_RHMaxSet.input)

        connect(itp_occ_ratio.output, occ_ratio)
        connect(itp_occ_active.output, occ_active)
        connect(itp_people_count.output, people_count)
        connect(itp_vent_min.output, vent_min)
        connect(itp_QOccSen.output, QOccSen)
        connect(itp_QOccLat.output, QOccLat)
        connect(itp_PLight.output, PLight)
        connect(itp_QLightRad.output, QLightRad)
        connect(itp_QLightCon.output, QLightCon)
        connect(itp_THeaSet.output, THeaSet)
        connect(itp_TCooSet.output, TCooSet)
        connect(itp_RHMinSet.output, RHMinSet)
        connect(itp_RHMaxSet.output, RHMaxSet)
    ]

    return ODESystem(eqs, t; name,
        systems = [
            clk,
            itp_occ_ratio, itp_occ_active, itp_people_count, itp_vent_min,
            itp_QOccSen, itp_QOccLat, itp_PLight, itp_QLightRad, itp_QLightCon,
            itp_THeaSet, itp_TCooSet, itp_RHMinSet, itp_RHMaxSet,
            occ_ratio, occ_active, people_count, vent_min,
            QOccSen, QOccLat, PLight, QLightRad, QLightCon,
            THeaSet, TCooSet, RHMinSet, RHMaxSet,
        ])
end


function ScheduleBus(; name::Symbol = :ScheduleBus,
                     time_col::Symbol = :time,
                     df::Union{Nothing,AbstractDataFrame} = nothing,
                     interp_method = ConstantInterpolation,
                     kwargs...)
    df_ = df

    if df_ === nothing
        for (_, v) in kwargs
            if v isa AbstractDataFrame
                df_ = v
                break
            end
        end
    end

    df_ === nothing && error("ScheduleBus: missing DataFrame. Pass `ScheduleBus(df)` or `ScheduleBus(df=...)`.")
    return ScheduleBus(DataFrame(df_); name = name, time_col = time_col, interp_method = interp_method)
end
