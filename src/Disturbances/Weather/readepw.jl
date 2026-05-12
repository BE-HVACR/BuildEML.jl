using Dates, CSV, DataFrames

function _epw_time_seconds(df_raw::DataFrame)
    years = Int.(df_raw.year)
    months = Int.(df_raw.month)
    days = Int.(df_raw.day)
    hours = Int.(df_raw.hour)
    minutes = Int.(df_raw.minute)
    is_hourly = all(m -> (m == 0 || m == 60), minutes)

    n = nrow(df_raw)
    times = Vector{Float64}(undef, n)
    base_year = years[1]
    t0 = DateTime(base_year, 1, 1, 0, 0, 0)
    nominal_year = 365.0 * 24.0 * 3600.0
    prev_t = -Inf
    cycle_shift = 0.0

    for i in 1:n
        y = years[i]
        mon = months[i]
        day = days[i]
        hour = hours[i]
        minute = minutes[i]

        if is_hourly
            dt = DateTime(y, mon, day, 0, 0, 0) + Hour(hour)
        else
            if minute == 60
                dt = DateTime(y, mon, day, 0, 0, 0) + Hour(hour)
            else
                dt = DateTime(y, mon, day, 0, 0, 0) + Hour(hour - 1) + Minute(minute)
            end
        end

        ti = Dates.value(dt - t0) / 1000 + cycle_shift
        while i > 1 && ti < prev_t
            cycle_shift += nominal_year
            ti += nominal_year
        end

        times[i] = ti
        prev_t = ti
    end

    return times, is_hourly
end

function _apply_reader_tmy3_hourly_alignment(df::DataFrame)
    # MBL ReaderTMY3 convention: h_k data is placed at t = k × 3600 (end-of-hour),
    # with the first record duplicated at t=0 as a periodic anchor.
    first_row = deepcopy(df[1, :])
    last_row  = deepcopy(df[end, :])
    df_shifted = vcat(DataFrame(first_row), df[1:end-1, :])
    df_shifted.time = collect(0.0:3600.0:(3600.0 * (nrow(df_shifted) - 1)))

    # Use the last EPW row (Dec 31 h24) as the periodic endpoint so that
    # interpolation at t = 8760×3600 stays near year-end values.
    endpoint = last_row
    endpoint.time = df_shifted.time[end] + 3600.0
    push!(df_shifted, endpoint)

    return df_shifted
end

"""
    ReadEPW(epw_path; return_raw=false)

Read an EPW file and return a processed DataFrame with standardized names/units.

Returned columns (units in brackets):
- time      [s]    : EPW-decoded local clock time in seconds.
                    For hourly EPW files, this function emulates
                    ReaderTMY3 alignment at year boundaries:
                    duplicate first record at t=0, drop year-end midnight record,
                    then append a periodic endpoint at t=end+dt.
- TDryBul   [K]    : dry-bulb temperature
- TDewPoi   [K]    : dew-point temperature
- relHum    [-]    : relative humidity (0-1)
- HumRat    [kg/kg]: humidity ratio
- pAtm      [Pa]   : atmospheric pressure
- TWetBul   [K]    : wet-bulb temperature (Stull 2011 approximation)
- HGloHor   [W/m2] : global horizontal irradiance (EPW :ghi)
- HDifHor   [W/m2] : diffuse horizontal irradiance (EPW :dhi)
- HDirNor   [W/m2] : direct normal irradiance (EPW :dni)
- HHorIR    [W/m2] : horizontal infrared irradiation (EPW :ghi_infrared)
- albedo    [-]    : ground reflectance (EPW :albedo)
- winDir    [rad]  : wind direction (deg -> rad)
- winSpe    [m/s]  : wind speed
"""
function ReadEPW(epw_path::AbstractString; return_raw::Bool = false)
    colnames = [:year, :month, :day, :hour, :minute, :data_source_unct,
        :temp_air, :temp_dew, :relative_humidity,
        :atmospheric_pressure, :etr, :etrn, :ghi_infrared, :ghi,
        :dni, :dhi, :global_hor_illum, :direct_normal_illum,
        :diffuse_horizontal_illum, :zenith_luminance,
        :wind_direction, :wind_speed, :total_sky_cover,
        :opaque_sky_cover, :visibility, :ceiling_height,
        :present_weather_observation, :present_weather_codes,
        :precipitable_water, :aerosol_optical_depth, :snow_depth,
        :days_since_last_snowfall, :albedo,
        :liquid_precipitation_depth, :liquid_precipitation_quantity]

    df_raw = DataFrame(CSV.File(epw_path; header = colnames, skipto = 9))
    return_raw && return df_raw

    time, is_hourly = _epw_time_seconds(df_raw)

    clamp01(x) = x < 0 ? 0.0 : (x > 1 ? 1.0 : x)
    deg2rad(x) = x * (pi / 180.0)

    TDryBul = Float64.(df_raw.temp_air) .+ 273.15
    TDewPoi = Float64.(df_raw.temp_dew) .+ 273.15
    relHum = clamp01.(Float64.(df_raw.relative_humidity) ./ 100.0)
    pAtm = Float64.(df_raw.atmospheric_pressure)

    TWetBul = wetbulb_TRH.(TDryBul, relHum)
    HumRat = w_TRHp.(TDryBul, relHum, pAtm)

    HGloHor = Float64.(df_raw.ghi)
    HDifHor = Float64.(df_raw.dhi)
    HDirNor = Float64.(df_raw.dni)
    HHorIR = Float64.(df_raw.ghi_infrared)
    albedo = Float64.(df_raw.albedo)
    winDir = deg2rad.(Float64.(df_raw.wind_direction))
    winSpe = Float64.(df_raw.wind_speed)

    df = DataFrame(;
        time, TDryBul, TDewPoi, relHum, pAtm, TWetBul, HumRat,
        HGloHor, HDifHor, HDirNor, HHorIR, albedo, winDir, winSpe,
    )

    if is_hourly
        df = _apply_reader_tmy3_hourly_alignment(df)
    else
        dt = nrow(df) > 1 ? (df.time[end] - df.time[end - 1]) : 3600.0
        last = deepcopy(df[end, :])
        last.time = df.time[end] + dt
        push!(df, last)
    end

    return df
end

const WEATHER_UNITS = Dict(
    :time => "s",
    :TDryBul => "K",
    :TDewPoi => "K",
    :relHum => "-",
    :pAtm => "Pa",
    :TWetBul => "K",
    :HumRat => "kg/kg_da",
    :HGloHor => "W/m^2",
    :HDifHor => "W/m^2",
    :HDirNor => "W/m^2",
    :HHorIR => "W/m^2",
    :albedo => "-",
    :winDir => "rad",
    :winSpe => "m/s",
)
