"""
    AkimaSpline(u, t)

Akima spline interpolation with linear extrapolation outside `t`.
The degenerate case is taken on an exact zero denominator, as the Akima method
specifies, rather than on a tolerance.

Usage: `itp = AkimaSpline(u, t); itp(x)`.
"""
struct AkimaSpline
    t::Vector{Float64}
    u::Vector{Float64}
    c0::Vector{Float64}   # per-interval cubic coefficients
    c1::Vector{Float64}
    c2::Vector{Float64}
end

function AkimaSpline(u::AbstractVector, t::AbstractVector)
    n = length(t)
    n == length(u) || error("AkimaSpline: length(u) != length(t).")
    n >= 3 || error("AkimaSpline: at least 3 points required.")
    all(diff(t) .> 0) || error("AkimaSpline: t must be strictly increasing.")

    tt = Float64.(t)
    uu = Float64.(u)

    # interval slopes, with two extra entries at each end
    dd = Vector{Float64}(undef, n + 3)
    for k in 1:(n - 1)
        dd[k + 2] = (uu[k + 1] - uu[k]) / (tt[k + 1] - tt[k])
    end

    # ghost slopes
    dd[1] = 3 * dd[3] - 2 * dd[4]
    dd[2] = 2 * dd[3] - dd[4]
    dd[n + 2] = 2 * dd[n + 1] - dd[n]
    dd[n + 3] = 3 * dd[n + 1] - 2 * dd[n]

    # left boundary slope
    c2cur = abs(dd[4] - dd[3]) + abs(dd[2] - dd[1])
    if c2cur > 0
        a = abs(dd[2] - dd[1]) / c2cur
        c2cur = (1 - a) * dd[2] + a * dd[3]
    else
        c2cur = 0.5 * dd[2] + 0.5 * dd[3]
    end

    c0 = Vector{Float64}(undef, n - 1)
    c1 = Vector{Float64}(undef, n - 1)
    c2 = Vector{Float64}(undef, n - 1)
    for k in 1:(n - 1)
        dx = tt[k + 1] - tt[k]
        c2[k] = c2cur
        c2cur = abs(dd[k + 4] - dd[k + 3]) + abs(dd[k + 2] - dd[k + 1])
        if c2cur > 0
            a = abs(dd[k + 2] - dd[k + 1]) / c2cur
            c2cur = (1 - a) * dd[k + 2] + a * dd[k + 3]
        else
            c2cur = 0.5 * dd[k + 2] + 0.5 * dd[k + 3]
        end
        c1[k] = (3 * dd[k + 2] - 2 * c2[k] - c2cur) / dx
        c0[k] = (c2[k] + c2cur - 2 * dd[k + 2]) / (dx * dx)
    end

    return AkimaSpline(tt, uu, c0, c1, c2)
end

function (itp::AkimaSpline)(x)
    t = itp.t
    n = length(t)
    if x < t[1]
        # linear extrapolation, left
        return itp.u[1] + itp.c2[1] * (x - t[1])
    elseif x >= t[n]
        # linear extrapolation, right
        v = t[n] - t[n - 1]
        der = (3 * itp.c0[n - 1] * v + 2 * itp.c1[n - 1]) * v + itp.c2[n - 1]
        return itp.u[n] + der * (x - t[n])
    else
        k = min(searchsortedlast(t, x), n - 1)
        v = x - t[k]
        return itp.u[k] + ((itp.c0[k] * v + itp.c1[k]) * v + itp.c2[k]) * v
    end
end

"""
    param_linear_interp(t, y_values, t_values)

Evaluate the linearly interpolated value at `t` from samples `y_values` defined on the grid `t_values`.
"""
function param_linear_interp(t, y_values, t_values)
    itp = LinearInterpolation(y_values, t_values; extrapolation = ExtrapolationType.Linear)
    return itp(t)
end

"""
    param_constant_interp(t, y_values, t_values)

Evaluate the piecewise-constant interpolated value at `t` from samples `y_values` defined on the grid `t_values`.
"""
function param_constant_interp(t, y_values, t_values)
    itp = ConstantInterpolation(y_values, t_values; extrapolation = ExtrapolationType.Constant)
    return itp(t)
end

@register_symbolic param_linear_interp(t, y::AbstractVector, t_grid::AbstractVector)
@register_symbolic param_constant_interp(t, y::AbstractVector, t_grid::AbstractVector)

"""
    ParameterizedSource(t_grid_input; if_constant_interpolation=false, name)

Create a source component whose output is driven by parameterized values sampled on the fixed time grid `t_grid_input`.

Set `if_constant_interpolation=true` to use piecewise-constant interpolation; otherwise the source uses linear interpolation.
"""
function ParameterizedSource(t_grid_input::AbstractVector; if_constant_interpolation::Bool = false, name)
    len = length(t_grid_input)

    @parameters y_params[1:len]
    @parameters t_grid_params[1:len]

    @named clk = ContinuousClock()
    @named output = RealOutput()
    if if_constant_interpolation
        eqs = [output.u ~ param_constant_interp(clk.output.u, y_params, t_grid_params)]
    else
        eqs = [output.u ~ param_linear_interp(clk.output.u, y_params, t_grid_params)]
    end

    defs = Dict{Any, Any}()
    for (i, val) in enumerate(t_grid_input)
        defs[t_grid_params[i]] = val
        defs[y_params[i]] = 0.0
    end

    ODESystem(eqs, t, [], [y_params..., t_grid_params...];
        systems = [clk, output],
        defaults = defs,
        name = name)
end

"""
    FirstOrderLag

A first-order lag model with time constant `tau`, input `u(t)`, and output `y(t)`.
"""
@mtkmodel FirstOrderLag begin
    @parameters begin
        tau = 30.0
    end
    @variables begin
        u(t)
        y(t)
    end
    @equations begin
        D(y) ~ (u - y) / tau
    end
end
