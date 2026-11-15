module Media

using ModelingToolkit
using ModelingToolkit: t_nounits as t

export WaterPort
include("Water/waterport.jl")
export cp_water, rho_water
include("Water/water_functions.jl")

export AirPort
include("Air/airport.jl")
export cp_da, cp_vapor, h_fg, k_mair, R_da, rho_air_ref, p_air_ref, patm
export h_g, h_da, h_f, h_Tw, rho_Tw, rho_Twp, rho_p, T_hw, p_saturation, w_TRHp, wetbulb_TRH, cp_moistair_w, w_sat, T_dew
include("Air/moistair_functions.jl")

end
