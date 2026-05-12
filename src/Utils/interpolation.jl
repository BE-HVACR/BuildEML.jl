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
