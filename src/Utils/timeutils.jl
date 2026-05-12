using Dates

"""
    dt_seconds(month, day; hour=0.0, base_year=2001)

Convert `(month, day, hour)` into seconds since `January 1, 00:00` of
`base_year`. A non-leap base year is preferred when a fixed 365-day calendar
is desired.
"""
function dt_seconds(month::Int, day::Int; hour::Real = 0.0, base_year::Int = 2001)
    h = floor(Int, hour)
    m = floor(Int, (hour - h) * 60)
    s = round(Int, ((hour - h) * 60 - m) * 60)
    t0 = DateTime(base_year, 1, 1, 0, 0, 0)
    dt = DateTime(base_year, month, day, h, m, s)
    return Dates.value(dt - t0) / 1000
end

"""
    build_tspan(start_month, start_day, end_month, end_day;
                start_hour=0.0, end_hour=24.0, inclusive_end=true, base_year=2001)

Return `(t0, t1)` in seconds since `January 1, 00:00` of `base_year`.
A non-leap base year is preferred when a fixed 365-day calendar is desired.

Rules for `inclusive_end`:
- `true`: `t1` is exactly the specified end timestamp.
- `false`: if `end_hour` is an integer hour, `t1` advances to the next hour
  mark; otherwise `t1` stays at the specified fractional end time.
"""
function build_tspan(start_month::Int, start_day::Int, end_month::Int, end_day::Int;
    start_hour::Real = 0.0, end_hour::Real = 24.0,
    inclusive_end::Bool = true, base_year::Int = 2001)

    t0 = dt_seconds(start_month, start_day; hour = start_hour, base_year = base_year)
    t1 = dt_seconds(end_month, end_day; hour = end_hour, base_year = base_year)

    if !inclusive_end && isapprox(end_hour, round(end_hour); atol = 1e-12)
        t1 += 3600
    end

    if t1 < t0
        throw(ArgumentError("Hi! BuildEML hasn't unlocked time-travel features yet: end time (t1=$t1) precedes start time (t0=$t0)."))
    end
    return (t0, t1)
end

#=
# Example checks
build_tspan(4, 12, 8, 20)
build_tspan(1, 1, 12, 31)[2] == 8759 * 3600
build_tspan(1, 1, 1, 1; start_hour = 0, end_hour = 1)
build_tspan(1, 1, 1, 1; start_hour = 0, end_hour = 1.5)
build_tspan(1, 1, 2, 30) # error
build_tspan(1, 1, 1, 1; start_hour = 1, end_hour = 0) # error
=#
