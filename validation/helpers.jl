using DataInterpolations

"""
    read_timetable(src, name)

Time and value columns of the `Modelica.Blocks.Sources.TimeTable` named `name`
in the Modelica source text `src`.
"""
function read_timetable(src, name)
    i = findfirst("TimeTable " * name * "(", src)
    i === nothing && error("table $name not found")
    k = findnext("table=[", src, last(i))
    j = findnext("])", src, last(k))
    rows = split(replace(src[(last(k) + 1):(first(j) - 1)], r"\s" => ""), ";")
    tv = [parse.(Float64, split(r, ",")) for r in rows]
    return first.(tv), last.(tv)
end

"""
    ramped(ts, vs; dt = 1e-3)

Linear interpolation of a TimeTable whose steps are written as repeated time
points, with each step turned into a linear ramp of `dt` seconds. Returns the
interpolation and its time points.
"""
function ramped(ts, vs; dt = 1e-3)
    tt = copy(ts)
    for i in 2:length(tt)
        tt[i] <= tt[i - 1] && (tt[i] = tt[i - 1] + dt)
    end
    return LinearInterpolation(vs, tt; extrapolation = ExtrapolationType.Constant), tt
end

"""
    read_mbl_reference(path)

Time grid and variables of a BuildingsPy reference result file
(MBL Resources/ReferenceResults/Dymola). The file stores `time=[t_start, t_end]`
and each variable on equidistant points over that span; a constant variable
is stored as two values and is expanded to the full length.
"""
function read_mbl_reference(path)
    ref = Dict{String, Vector{Float64}}()
    for l in readlines(path)
        k = findfirst("=[", l)
        (k === nothing || !endswith(l, "]")) && continue
        ref[l[1:(first(k) - 1)]] = parse.(Float64, strip.(split(l[(last(k) + 1):(end - 1)], ",")))
    end
    t_start, t_end = ref["time"]
    n = maximum(length(v) for (name, v) in ref if name != "time")
    for (name, v) in ref
        name != "time" && length(v) == 2 && (ref[name] = fill(v[1], n))
    end
    return collect(range(t_start, t_end, length = n)), ref
end
