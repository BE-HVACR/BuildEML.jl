function _expand_hourly_profile(hourly_values::AbstractVector, bin_seconds::Real, label::AbstractString)
    length(hourly_values) == 24 || throw(ArgumentError("$label hourly profile must contain exactly 24 values."))
    bin_seconds > 0 || throw(ArgumentError("bin_seconds must be positive."))

    bins_per_hour = 3600.0 / float(bin_seconds)
    isapprox(bins_per_hour, round(bins_per_hour); atol = 1e-9, rtol = 1e-9) ||
        throw(ArgumentError("MakeScheduleTypicalOffice requires bin_seconds to evenly divide one hour. Got bin_seconds=$(bin_seconds)."))

    repeat_count = Int(round(bins_per_hour))
    return repeat(Float64.(hourly_values), inner = repeat_count)
end


"""
    MakeScheduleTypicalOffice(; ...)

Build a typical office operation schedule using the ASHRAE 90.1-1989 office
occupancy profile:
- Weekday equivalent full-load hours: 9.2
- Saturday equivalent full-load hours: 2.0
- Sunday/holiday equivalent full-load hours: 0.6

This helper delegates to `MakeSchedule(...)` with the following reference-based
defaults. These are intended as estimated typical-office assumptions rather than
one-to-one mandatory values for every office:
- weekday / saturday / sunday occupancy profiles pre-filled
- `occupant_density = 0.05` [person/m2]
  from ASHRAE 62.1-2016 Table 6.2.2.1, office buildings - occupant density
- `min_vent_rate_perperson = 0.0085` [m3/s-person]
  converted from 8.5 L/s-person, from ASHRAE 62.1-2016 Table 6.2.2.1,
  office buildings - office space Combined Outdoor Air Rate
- `occ_sen_heatgain_perperson = 75.0` [W/person]
  reference occupant sensible heat gain
- `occ_lat_heatgain_perperson = 45.0` [W/person]
  reference occupant latent heat gain, referenced from: ASHRAE 55-2017 Thermal Environmental Conditions for Human Occupancy, 2017.
- `if_other_loads = true`
- lighting ratios follow the same weekday / saturday / sunday profiles as occupancy
- `lighting_density = 7.5347` [W/m2]
  converted from 0.7 W/ft2, referenced from:
  Myer M, Seeger K. Lighting Changes in ASHRAE/IES Standard 90.1-2022. ASHRAE Journal 2023;65.
- `lighting_radiant_ratio = 0.4`
  reference estimate for recessed LED troffers, referenced from:
  Liu R, Zhou X, Lochhead SJ, Zhong Z, Huynh CV, Maxwell GM.
  Low-energy LED lighting heat gain distribution in buildings, part II:
  LED luminaire selection and test results.
  Science and Technology for the Built Environment 2017;23:688-708.
- `thea_set_occ = 21 degC`, `thea_set_uno = 18 degC`,
  `tcoo_set_occ = 24 degC`, `tcoo_set_uno = 27 degC`
  reference setback temperatures from:
  Elehwany H, Gunay B, Ouf M, Cotrufo N, Venne J-S.
  Evaluating common supply air temperature setpoint reset strategies with varying occupancy patterns and behaviours.
  Building and Environment 2024;266:112129.
- `rhmin_set_occ = 0.40`, `rhmin_set_uno = 0.20`,
  `rhmax_set_occ = 0.60`, `rhmax_set_uno = 0.70`
  reference humidity limits from:
  Zhivov A, Rose W, Patenaude R, Williams WJ.
  Requirements for Building Thermal Conditions under Normal and Emergency Operations in Extreme Climates.
  ASHRAE Journal. 2021 Jan;127(1):693-704.
  The unoccupied lower-bound default of 20% is an engineering assumption layered on top of that reference.

You still must provide the remaining required operation inputs, such as:
- `n_people_full`, or `floor_area`
"""
function MakeScheduleTypicalOffice(;
        n_people_full::Union{Nothing,Real} = nothing,
        floor_area::Real,
        occupant_density::Real = 0.05,
        min_vent_rate_perperson::Real = 0.0085,
        occ_sen_heatgain_perperson::Real = 75.0,
        occ_lat_heatgain_perperson::Real = 45.0,
        if_other_loads::Bool = true,
        lighting_density::Real = 0.7 / 0.09290304,
        lighting_radiant_ratio::Real = 0.4,
        thea_set_occ::Real = unit_C2K(21.0),
        thea_set_uno::Real = unit_C2K(18.0),
        tcoo_set_occ::Real = unit_C2K(24.0),
        tcoo_set_uno::Real = unit_C2K(27.0),
        rhmin_set_occ::Real = 0.40,
        rhmin_set_uno::Real = 0.20,
        rhmax_set_occ::Real = 0.60,
        rhmax_set_uno::Real = 0.70,
        if_us_holidays::Bool = false,
        explicit_day_types = nothing,
        if_occupant_heatgain::Bool = true,
        bin_seconds::Real = 3600.0,
        base_year::Int = 2001,
        occ_eps::Real = 1e-6)
    occ_ratio_weekday = _expand_hourly_profile(
        [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.1, 0.2,
            0.95, 0.95, 0.95, 0.95,
            0.5,
            0.95, 0.95, 0.95, 0.95,
            0.3,
            0.1, 0.1, 0.1, 0.1,
            0.05, 0.05,
        ],
        bin_seconds,
        "weekday",
    )
    occ_ratio_saturday = _expand_hourly_profile(
        [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.1, 0.1,
            0.3, 0.3, 0.3, 0.3,
            0.1, 0.1, 0.1, 0.1, 0.1,
            0.05, 0.05,
            0.0, 0.0, 0.0, 0.0, 0.0,
        ],
        bin_seconds,
        "saturday",
    )
    occ_ratio_sunday = _expand_hourly_profile(
        [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ],
        bin_seconds,
        "sunday/holiday",
    )

    return MakeSchedule(;
        occ_ratio_weekday = occ_ratio_weekday,
        occ_ratio_saturday = occ_ratio_saturday,
        occ_ratio_sunday = occ_ratio_sunday,
        explicit_day_types = explicit_day_types,
        if_us_holidays = if_us_holidays,
        n_people_full = n_people_full,
        floor_area = floor_area,
        occupant_density = occupant_density,
        min_vent_rate_perperson = min_vent_rate_perperson,
        if_occupant_heatgain = if_occupant_heatgain,
        occ_sen_heatgain_perperson = occ_sen_heatgain_perperson,
        occ_lat_heatgain_perperson = occ_lat_heatgain_perperson,
        if_other_loads = if_other_loads,
        lighting_ratio_weekday = occ_ratio_weekday,
        lighting_ratio_saturday = occ_ratio_saturday,
        lighting_ratio_sunday = occ_ratio_sunday,
        lighting_density = lighting_density,
        lighting_radiant_ratio = lighting_radiant_ratio,
        thea_set_occ = thea_set_occ,
        thea_set_uno = thea_set_uno,
        tcoo_set_occ = tcoo_set_occ,
        tcoo_set_uno = tcoo_set_uno,
        rhmin_set_occ = rhmin_set_occ,
        rhmin_set_uno = rhmin_set_uno,
        rhmax_set_occ = rhmax_set_occ,
        rhmax_set_uno = rhmax_set_uno,
        occ_eps = occ_eps,
        bin_seconds = bin_seconds,
        base_year = base_year,
    )
end
