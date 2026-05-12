module Utils

using Dates
using Unitful
using ModelingToolkit, Symbolics, IfElse
using ModelingToolkitStandardLibrary.Blocks
using DataInterpolations
using ModelingToolkit: t_nounits as t, D_nounits as D

export smooth_H, smooth_xH, smooth_max, smooth_min, smooth_abs, smooth_sign, smooth_clamp, smooth_clamp01, smooth_on_eps, snapconst, α_soft, soft_blend
include("smooth.jl")

export build_tspan, dt_seconds
include("timeutils.jl")

export unit_K2C, unit_C2K, to_SI, to_SI_T, to_SI_p, to_SI_massflow, to_SI_volflow, to_SI_energy, to_SI_power, ensure_SI
include("units.jl")

export ParameterizedSource, FirstOrderLag
include("interpolation.jl")

end
