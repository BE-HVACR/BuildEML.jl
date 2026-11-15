# Moist air properties on a dry-air mass basis.
#
# Refer to the MBL v11:
# - p_air_ref, rho_air_ref, rho_p: Buildings.Media.Air reference_p, dStp, and its density law
# - cp_da, cp_vapor, h_fg, k_mair: Buildings.Utilities.Psychrometrics.Constants
# - R_da: MSL v4 Modelica.Media.IdealGases.Common.SingleGasesData.Air R_s
# - h_g, h_da, h_f, h_Tw, cp_moistair_w, T_hw: Buildings.Media.Air enthalpy relations
# - w_TRHp, w_sat: Buildings.Utilities.Psychrometrics.Functions.X_pSatpphi
# - p_saturation, T_dew: Buildings.Utilities.Psychrometrics.Functions.saturationPressureLiquid,
#   of which T_dew is the inverse
# - wetbulb_TRH: Buildings.Utilities.Psychrometrics.TWetBul_TDryBulPhi
#   with approximateWetBulb = true
# - wetbulb_balance_Twp: Buildings.Utilities.Psychrometrics.TWetBul_TDryBulXi
#   with approximateWetBulb = false (wetbulb_Twp solves it with Newton iterations)
# Not from MBL:
# - rho_Twp

const p_air_ref = 101325.0   # Pa
const rho_air_ref = 1.2      # kg/m3
const patm = 101325.0        # Pa, standard atmosphere
const cp_da = 1006.0         # J/kg_da/K
const cp_vapor = 1860.0      # J/kg_w/K
const h_fg = 2501014.5       # J/kg_w, latent heat of vaporization at 0 degC
const k_mair = 0.6219647130774989   # [-] molar mass ratio, water vapor / dry air
const R_da = 287.0512249529787      # J/kg_da/K

# enthalpy
h_g(T) = h_fg + cp_vapor * (T - 273.15)   # J/kg_w
h_da(T) = cp_da * (T - 273.15)            # J/kg_da
h_f(T) = cp_water * (T - 273.15)          # J/kg_w

h_Tw(T, w) = cp_da * (T - 273.15) + w * (h_fg + cp_vapor * (T - 273.15))  # J/kg_da
cp_moistair_w(w) = cp_da + w * cp_vapor  # dh/dT on a dry-air basis
T_hw(h, w) = (h - h_fg * w) / (cp_da + cp_vapor * w) + 273.15  # K

# density
rho_Twp(T, w, p) = p / (R_da * T * (1 + w / k_mair))  # kg_da/m3
rho_p(p) = p * rho_air_ref / p_air_ref

# humidity and saturation
w_TRHp(T, RH, p) = k_mair * (RH * p_saturation(T)) / (p - (RH * p_saturation(T)))  # kg_w/kg_da

p_saturation(T) = 611.657 * exp(17.2799 - 4102.99 / (T - 35.719))  # Pa

wetbulb_TRH(T, RH) = 273.15 + (T - 273.15) * atan(0.151977 * sqrt(100 * RH + 8.313659)) +
    atan((T - 273.15) + 100 * RH) - atan(100 * RH - 1.676331) +
    0.00391838 * (100 * RH)^(3 / 2) * atan(0.023101 * 100 * RH) - 4.686035  # K

w_sat(T) = k_mair * p_saturation(T) / (patm - p_saturation(T))  # kg_w/kg_da

T_dew(w) = let pw = w * p_air_ref / (k_mair + w)
    35.719 + 4102.99 / (17.2799 - log(pw / 611.657))
end  # K

# wetbulb_balance_Twp(T_wb, T, w, p) = 0 defines the wet-bulb temperature implicitly
function wetbulb_balance_Twp(T_wb, T, w, p)
    psat = p_saturation(T_wb)
    w_sb = k_mair * psat / (p - psat)
    return (T_wb - 273.15) * (cp_da + w_sb * cp_vapor + (w - w_sb) * cp_water) -
           ((T - 273.15) * (cp_da + w * cp_vapor) + (w - w_sb) * h_fg)
end

function wetbulb_Twp(T, w, p = patm)
    T_wb = T
    for _ in 1:50
        r = wetbulb_balance_Twp(T_wb, T, w, p)
        dr = (wetbulb_balance_Twp(T_wb + 1e-6, T, w, p) - r) / 1e-6
        step = r / dr
        T_wb -= step
        abs(step) < 1e-10 && return T_wb
    end
    error("wetbulb_Twp did not converge for T=$(T), w=$(w), p=$(p)")
end  # K
