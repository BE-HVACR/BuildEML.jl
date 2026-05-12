module Media

using ModelingToolkit
using ModelingToolkit: t_nounits as t

export AirPort
include("Air/airport.jl")
export cp_da
export h_g, h_da, h_f, h_Tw, rho_Tw, rho_Twp, T_hw, patm, p_saturation, w_TRHp, wetbulb_TRH, cp_moistair_w
include("Air/moistair_functions.jl")

export WaterPort
include("Water/waterport.jl")
export cp_water, rho_water
include("Water/water_functions.jl")

end
