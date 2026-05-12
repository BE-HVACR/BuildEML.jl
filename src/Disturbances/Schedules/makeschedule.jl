const SCHEDULE_REQUIRED_COLUMNS = (
    :occRatio,
    :occActive,
    :peopleCount,
    :ventMin,
    :PLight,
    :QLightRad,
    :QLightCon,
    :THeaSet,
    :TCooSet,
    :RHMinSet,
    :RHMaxSet,
    :QOccSen,
    :QOccLat,
)


"""
    weekly_schedule_values(weekday; saturday=weekday, sunday=saturday)

Expand three day profiles into one weekly profile:
- 5 x weekday
- 1 x saturday
- 1 x sunday
"""
function weekly_schedule_values(weekday::AbstractVector;
                                saturday::AbstractVector = weekday,
                                sunday::AbstractVector = saturday)
    n = length(weekday)
    n > 0 || throw(ArgumentError("weekday profile must contain at least one entry."))
    length(saturday) == n || throw(ArgumentError("saturday profile must have the same length as weekday profile."))
    length(sunday) == n || throw(ArgumentError("sunday profile must have the same length as weekday profile."))
    return vcat(weekday, weekday, weekday, weekday, weekday, saturday, sunday)
end


function _validate_day_profile(profile::AbstractVector, bin_seconds::Real, label::AbstractString)
    bin_seconds > 0 || throw(ArgumentError("bin_seconds must be positive."))
    length(profile) > 0 || throw(ArgumentError("$label profile must contain at least one entry."))

    day_seconds = length(profile) * float(bin_seconds)
    isapprox(day_seconds, 24.0 * 3600.0; atol = 1e-9, rtol = 1e-9) ||
        throw(ArgumentError("$label profile must cover exactly 24 hours. " *
                            "Got $(length(profile)) bins with bin_seconds=$(bin_seconds)."))
end


function _default_unoccupied_profile_like(profile::AbstractVector)
    return zeros(Float64, length(profile))
end


function _observed_fixed_holiday(year::Int, month::Int, day::Int)
    d = Date(year, month, day)
    dow = dayofweek(d)
    if dow == 6
        return d - Day(1)
    elseif dow == 7
        return d + Day(1)
    else
        return d
    end
end


function _nth_weekday_of_month(year::Int, month::Int, weekday::Int, nth::Int)
    d = Date(year, month, 1)
    while dayofweek(d) != weekday
        d += Day(1)
    end
    return d + Day(7 * (nth - 1))
end


function _last_weekday_of_month(year::Int, month::Int, weekday::Int)
    d = Date(year, month, daysinmonth(Date(year, month, 1)))
    while dayofweek(d) != weekday
        d -= Day(1)
    end
    return d
end


# U.S. holidays are interpreted using federal holidays plus their observed dates.
function _us_federal_holidays_observed(year::Int)
    holidays = Date[
        _observed_fixed_holiday(year, 1, 1),   # New Year's Day
        _nth_weekday_of_month(year, 1, 1, 3),  # Martin Luther King Jr. Day
        _nth_weekday_of_month(year, 2, 1, 3),  # Washington's Birthday / Presidents Day
        _last_weekday_of_month(year, 5, 1),    # Memorial Day
        _nth_weekday_of_month(year, 9, 1, 1),  # Labor Day
        _nth_weekday_of_month(year, 10, 1, 2), # Columbus Day
        _observed_fixed_holiday(year, 11, 11), # Veterans Day
        _nth_weekday_of_month(year, 11, 4, 4), # Thanksgiving Day
        _observed_fixed_holiday(year, 7, 4),   # Independence Day
        _observed_fixed_holiday(year, 12, 25), # Christmas Day
    ]

    if year >= 2021
        push!(holidays, _observed_fixed_holiday(year, 6, 19)) # Juneteenth National Independence Day
    end

    return Set(holidays)
end


function _resolve_n_people_full(; n_people_full, floor_area, occupant_density)
    if n_people_full !== nothing
        return float(n_people_full)
    end

    occupant_density !== nothing ||
        throw(ArgumentError("Provide either `n_people_full`, or `occupant_density` together with required `floor_area`."))

    return float(floor_area) * float(occupant_density)
end


function _resolve_occupant_load_gain(; if_occupant_heatgain::Bool, gain, label::AbstractString)
    if !if_occupant_heatgain
        return 0.0
    end

    gain === nothing &&
        throw(ArgumentError("`$label` must be provided when `if_occupant_heatgain=true`."))

    return float(gain)
end


function _resolve_other_load_param(; if_other_loads::Bool, value, label::AbstractString)
    if !if_other_loads
        return 0.0
    end

    value === nothing &&
        throw(ArgumentError("`$label` must be provided when `if_other_loads=true`."))

    return float(value)
end


function _resolve_lighting_power_full(; if_other_loads::Bool, lighting_power_full, lighting_density, floor_area::Real)
    if !if_other_loads
        return 0.0
    end

    if lighting_power_full !== nothing
        return float(lighting_power_full)
    end

    lighting_density !== nothing ||
        throw(ArgumentError("Provide either `lighting_power_full`, or `lighting_density` together with required `floor_area` when `if_other_loads=true`."))

    return float(lighting_density) * float(floor_area)
end


function _build_schedule_rows(day_profile::AbstractVector,
                              lighting_profile::AbstractVector,
                              n_people_full::Real,
                              min_vent_rate_perperson::Real,
                              occ_sen_heatgain_perperson::Real,
                              occ_lat_heatgain_perperson::Real,
                              lighting_power_full::Real,
                              lighting_radiant_ratio::Real,
                              thea_set_occ::Real,
                              thea_set_uno::Real,
                              tcoo_set_occ::Real,
                              tcoo_set_uno::Real,
                              rhmin_set_occ::Real,
                              rhmin_set_uno::Real,
                              rhmax_set_occ::Real,
                              rhmax_set_uno::Real,
                              occ_eps::Real)
    occ_ratio = clamp.(Float64.(day_profile), 0.0, 1.0)
    occ_active = Float64.(occ_ratio .> occ_eps)
    lighting_ratio = clamp.(Float64.(lighting_profile), 0.0, 1.0)

    people_count = float(n_people_full) .* occ_ratio
    vent_min = float(min_vent_rate_perperson) .* people_count
    QOccSen = float(occ_sen_heatgain_perperson) .* people_count
    QOccLat = float(occ_lat_heatgain_perperson) .* people_count
    PLight = float(lighting_power_full) .* lighting_ratio
    QLightRad = float(lighting_radiant_ratio) .* PLight
    QLightCon = (1.0 - float(lighting_radiant_ratio)) .* PLight

    THeaSet = ifelse.(occ_active .> 0.5, float(thea_set_occ), float(thea_set_uno))
    TCooSet = ifelse.(occ_active .> 0.5, float(tcoo_set_occ), float(tcoo_set_uno))
    RHMinSet = ifelse.(occ_active .> 0.5, float(rhmin_set_occ), float(rhmin_set_uno))
    RHMaxSet = ifelse.(occ_active .> 0.5, float(rhmax_set_occ), float(rhmax_set_uno))

    return (; occ_ratio, occ_active, people_count, vent_min, QOccSen, QOccLat, PLight, QLightRad, QLightCon, THeaSet, TCooSet, RHMinSet, RHMaxSet)
end


function _normalize_explicit_day_types(explicit_day_types, base_year::Int)
    explicit_day_types === nothing && return nothing

    normalized = Dict{Date,Symbol}()
    for (k, v) in pairs(explicit_day_types)
        if k isa Date
            day = Date(base_year, month(k), day(k))
        elseif k isa Tuple && length(k) == 2
            month_i, day_i = k
            (month_i isa Integer && day_i isa Integer) ||
                throw(ArgumentError("Tuple keys of explicit_day_types must be `(month::Int, day::Int)`. Got $(repr(k))."))
            day = Date(base_year, Int(month_i), Int(day_i))
        else
            throw(ArgumentError("Keys of explicit_day_types must be `Date` or `(month, day)`. Got $(typeof(k))."))
        end
        v_sym = v isa Symbol ? v : Symbol(v)
        v_sym in (:workday, :saturday, :sunday, :holiday) ||
            throw(ArgumentError("Values of explicit_day_types must be `:workday`, `:saturday`, `:sunday`, or `:holiday`. Got $(repr(v))."))
        normalized[day] = v_sym
    end
    return normalized
end


"""
    MakeSchedule(; ...)

Build a full-year schedule table for one calendar year.

Inputs:
- daily occupancy-ratio profiles for weekday / saturday / sunday
  Each profile must satisfy `length(profile) * bin_seconds == 24 * 3600`.
  Example: 24 values for `bin_seconds = 3600`, 48 values for `bin_seconds = 1800`.
  Example hourly profile:
  `occ_ratio_weekday = vcat(zeros(7), ones(12), zeros(5))`
  By default, saturday and sunday profiles are fully unoccupied.
  Optional `explicit_day_types = Dict((1, 3) => :workday, (1, 4) => :holiday)`
  overrides the normal calendar logic: all unspecified dates become `holiday`, and only
  the listed month/day dates use the requested day type. Allowed values are
  `:workday`, `:saturday`, `:sunday`, and `:holiday`. Keys may be `(month, day)`
  tuples or `Date`; in both cases only the month/day is used and it is mapped onto
  `base_year`. U.S. holidays, when enabled, use the sunday profile.
- design parameters for people count, minimum ventilation, temperature setpoints,
  and humidity setpoints
  `n_people_full` directly sets the full-occupancy people count and overrides
  `floor_area * occupant_density` when provided.
  `floor_area` is required. Provide either `n_people_full`, or `occupant_density`.
  `min_vent_rate_perperson` is the minimum ventilation rate per person [m3/s].
  `if_occupant_heatgain=true` adds occupant sensible/latent loads and requires
  `occ_sen_heatgain_perperson` and `occ_lat_heatgain_perperson` [W/person].
  `if_other_loads=true` adds lighting loads and requires
  `lighting_ratio_weekday`, `lighting_radiant_ratio` [-], and either
  `lighting_power_full` [W] or `lighting_density` [W/m^2].
  `lighting_ratio_saturday` and `lighting_ratio_sunday` default to fully unoccupied.

Outputs:
- `time`      [s]
- `occRatio`  [-]
- `occActive` [-] 0 or 1
- `peopleCount`
- `ventMin`   [m3/s]
- `QOccSen`   [W]
- `QOccLat`   [W]
- `PLight`    [W]
- `QLightRad` [W]
- `QLightCon` [W]
- `THeaSet`   [K]
- `TCooSet`   [K]
- `RHMinSet`  [-]
- `RHMaxSet`  [-]

`base_year` determines the weekday alignment and, when enabled, the U.S. holiday dates.
Any year can be used. For leap years, the schedule is still forced to 365 days / 8760 h
by truncating the final day (`12/31`).
"""
function MakeSchedule(;
        occ_ratio_weekday::AbstractVector,
        occ_ratio_saturday::Union{Nothing,AbstractVector} = nothing,
        occ_ratio_sunday::Union{Nothing,AbstractVector} = nothing,
        lighting_ratio_weekday::Union{Nothing,AbstractVector} = nothing,
        lighting_ratio_saturday::Union{Nothing,AbstractVector} = nothing,
        lighting_ratio_sunday::Union{Nothing,AbstractVector} = nothing,
        explicit_day_types = nothing,
        if_us_holidays::Bool = false,
        n_people_full::Union{Nothing,Real} = nothing,
        floor_area::Real,
        occupant_density::Union{Nothing,Real} = nothing,
        min_vent_rate_perperson::Real,
        if_occupant_heatgain::Bool = true,
        occ_sen_heatgain_perperson::Union{Nothing,Real} = nothing,
        occ_lat_heatgain_perperson::Union{Nothing,Real} = nothing,
        if_other_loads::Bool = false,
        lighting_power_full::Union{Nothing,Real} = nothing,
        lighting_density::Union{Nothing,Real} = nothing,
        lighting_radiant_ratio::Union{Nothing,Real} = nothing,
        thea_set_occ::Real,
        thea_set_uno::Real,
        tcoo_set_occ::Real,
        tcoo_set_uno::Real,
        rhmin_set_occ::Real,
        rhmin_set_uno::Real,
        rhmax_set_occ::Real,
        rhmax_set_uno::Real,
        occ_eps::Real = 1e-6,
        bin_seconds::Real = 3600.0,
        base_year::Int = 2001)
    occ_ratio_saturday === nothing && (occ_ratio_saturday = _default_unoccupied_profile_like(occ_ratio_weekday))
    occ_ratio_sunday === nothing && (occ_ratio_sunday = _default_unoccupied_profile_like(occ_ratio_weekday))
    if if_other_loads
        lighting_ratio_weekday === nothing &&
            throw(ArgumentError("`lighting_ratio_weekday` must be provided when `if_other_loads=true`."))
        lighting_ratio_saturday === nothing && (lighting_ratio_saturday = _default_unoccupied_profile_like(lighting_ratio_weekday))
        lighting_ratio_sunday === nothing && (lighting_ratio_sunday = _default_unoccupied_profile_like(lighting_ratio_weekday))
    else
        lighting_ratio_weekday === nothing && (lighting_ratio_weekday = _default_unoccupied_profile_like(occ_ratio_weekday))
        lighting_ratio_saturday === nothing && (lighting_ratio_saturday = _default_unoccupied_profile_like(lighting_ratio_weekday))
        lighting_ratio_sunday === nothing && (lighting_ratio_sunday = _default_unoccupied_profile_like(lighting_ratio_weekday))
    end
    explicit_day_types = _normalize_explicit_day_types(explicit_day_types, base_year)
    n_people_full_effective = _resolve_n_people_full(;
        n_people_full = n_people_full,
        floor_area = floor_area,
        occupant_density = occupant_density,
    )
    occ_sen_heatgain_perperson_effective = _resolve_occupant_load_gain(;
        if_occupant_heatgain = if_occupant_heatgain,
        gain = occ_sen_heatgain_perperson,
        label = "occ_sen_heatgain_perperson",
    )
    occ_lat_heatgain_perperson_effective = _resolve_occupant_load_gain(;
        if_occupant_heatgain = if_occupant_heatgain,
        gain = occ_lat_heatgain_perperson,
        label = "occ_lat_heatgain_perperson",
    )
    lighting_power_full_effective = _resolve_lighting_power_full(;
        if_other_loads = if_other_loads,
        lighting_power_full = lighting_power_full,
        lighting_density = lighting_density,
        floor_area = floor_area,
    )
    lighting_radiant_ratio_effective = _resolve_other_load_param(;
        if_other_loads = if_other_loads,
        value = lighting_radiant_ratio,
        label = "lighting_radiant_ratio",
    )

    _validate_day_profile(occ_ratio_weekday, bin_seconds, "weekday")
    _validate_day_profile(occ_ratio_saturday, bin_seconds, "saturday")
    _validate_day_profile(occ_ratio_sunday, bin_seconds, "sunday")
    _validate_day_profile(lighting_ratio_weekday, bin_seconds, "lighting weekday")
    _validate_day_profile(lighting_ratio_saturday, bin_seconds, "lighting saturday")
    _validate_day_profile(lighting_ratio_sunday, bin_seconds, "lighting sunday")

    length(occ_ratio_saturday) == length(occ_ratio_weekday) ||
        throw(ArgumentError("saturday profile must have the same number of bins as weekday profile."))
    length(occ_ratio_sunday) == length(occ_ratio_weekday) ||
        throw(ArgumentError("sunday profile must have the same number of bins as weekday profile."))
    length(lighting_ratio_weekday) == length(occ_ratio_weekday) ||
        throw(ArgumentError("lighting weekday profile must have the same number of bins as weekday profile."))
    length(lighting_ratio_saturday) == length(occ_ratio_weekday) ||
        throw(ArgumentError("lighting saturday profile must have the same number of bins as weekday profile."))
    length(lighting_ratio_sunday) == length(occ_ratio_weekday) ||
        throw(ArgumentError("lighting sunday profile must have the same number of bins as weekday profile."))

    times = Float64[]
    occ_ratio_all = Float64[]
    occ_active_all = Float64[]
    people_count_all = Float64[]
    vent_min_all = Float64[]
    QOccSen_all = Float64[]
    QOccLat_all = Float64[]
    PLight_all = Float64[]
    QLightRad_all = Float64[]
    QLightCon_all = Float64[]
    THeaSet_all = Float64[]
    TCooSet_all = Float64[]
    RHMinSet_all = Float64[]
    RHMaxSet_all = Float64[]

    bins_per_day = length(occ_ratio_weekday)
    sizehint!(times, 365 * bins_per_day)
    sizehint!(occ_ratio_all, 365 * bins_per_day)
    sizehint!(occ_active_all, 365 * bins_per_day)
    sizehint!(people_count_all, 365 * bins_per_day)
    sizehint!(vent_min_all, 365 * bins_per_day)
    sizehint!(QOccSen_all, 365 * bins_per_day)
    sizehint!(QOccLat_all, 365 * bins_per_day)
    sizehint!(PLight_all, 365 * bins_per_day)
    sizehint!(QLightRad_all, 365 * bins_per_day)
    sizehint!(QLightCon_all, 365 * bins_per_day)
    sizehint!(THeaSet_all, 365 * bins_per_day)
    sizehint!(TCooSet_all, 365 * bins_per_day)
    sizehint!(RHMinSet_all, 365 * bins_per_day)
    sizehint!(RHMaxSet_all, 365 * bins_per_day)

    t_now = 0.0
    start_date = Date(base_year, 1, 1)
    # Keep the generated schedule at 365 days even in leap years.
    end_date = Dates.daysinyear(base_year) == 366 ? Date(base_year, 12, 30) : Date(base_year, 12, 31)
    holidays = if_us_holidays ? _us_federal_holidays_observed(base_year) : Set{Date}()

    for day in start_date:Day(1):end_date
        if explicit_day_types !== nothing
            day_type = get(explicit_day_types, day, :holiday)
            if day_type === :workday
                day_profile = occ_ratio_weekday
                lighting_profile = lighting_ratio_weekday
            elseif day_type === :saturday
                day_profile = occ_ratio_saturday
                lighting_profile = lighting_ratio_saturday
            else
                day_profile = occ_ratio_sunday
                lighting_profile = lighting_ratio_sunday
            end
        else
            dow = dayofweek(day)
            if if_us_holidays && (day in holidays)
                day_profile = occ_ratio_sunday
                lighting_profile = lighting_ratio_sunday
            else
                day_profile = dow <= 5 ? occ_ratio_weekday : (dow == 6 ? occ_ratio_saturday : occ_ratio_sunday)
                lighting_profile = dow <= 5 ? lighting_ratio_weekday : (dow == 6 ? lighting_ratio_saturday : lighting_ratio_sunday)
            end
        end

        rows = _build_schedule_rows(
            day_profile,
            lighting_profile,
            n_people_full_effective,
            min_vent_rate_perperson,
            occ_sen_heatgain_perperson_effective,
            occ_lat_heatgain_perperson_effective,
            lighting_power_full_effective,
            lighting_radiant_ratio_effective,
            thea_set_occ,
            thea_set_uno,
            tcoo_set_occ,
            tcoo_set_uno,
            rhmin_set_occ,
            rhmin_set_uno,
            rhmax_set_occ,
            rhmax_set_uno,
            occ_eps,
        )

        for i in eachindex(day_profile)
            push!(times, t_now)
            push!(occ_ratio_all, rows.occ_ratio[i])
            push!(occ_active_all, rows.occ_active[i])
            push!(people_count_all, rows.people_count[i])
            push!(vent_min_all, rows.vent_min[i])
            push!(QOccSen_all, rows.QOccSen[i])
            push!(QOccLat_all, rows.QOccLat[i])
            push!(PLight_all, rows.PLight[i])
            push!(QLightRad_all, rows.QLightRad[i])
            push!(QLightCon_all, rows.QLightCon[i])
            push!(THeaSet_all, rows.THeaSet[i])
            push!(TCooSet_all, rows.TCooSet[i])
            push!(RHMinSet_all, rows.RHMinSet[i])
            push!(RHMaxSet_all, rows.RHMaxSet[i])
            t_now += float(bin_seconds)
        end
    end

    return DataFrame(
        time = times,
        occRatio = occ_ratio_all,
        occActive = occ_active_all,
        peopleCount = people_count_all,
        ventMin = vent_min_all,
        QOccSen = QOccSen_all,
        QOccLat = QOccLat_all,
        PLight = PLight_all,
        QLightRad = QLightRad_all,
        QLightCon = QLightCon_all,
        THeaSet = THeaSet_all,
        TCooSet = TCooSet_all,
        RHMinSet = RHMinSet_all,
        RHMaxSet = RHMaxSet_all,
    )
end
