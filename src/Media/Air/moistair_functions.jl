const p_air_ref = 101325.0   # Pa
const patm = 101325.0   # Pa
const rho_air_ref = 1.2   # kg/m3
const cp_da = 1005.45   # J/kg_da/K

# ASHRAE Handbook - Fundamentals, Chapter 1 "Psychrometrics"
h_g(T) = 1e3 * (2501 + 1.86 * (T - 273.15))   # J/kg_w
h_da(T) = cp_da * (T - 273.15)                # J/kg_da
h_f(T) = 4186 * (T - 273.15)                  # J/kg_w

h_Tw(T, w) = cp_da * (T - 273.15) + w * 1e3 * (2501 + 1.86 * (T - 273.15))  # J/kg_da
cp_moistair_w(w) = cp_da + w * 1.86 * 1e3  # dh/dT on a dry-air basis
T_hw(h, w) = (h / 1e3 - 2501 * w) / (1.005 + 1.86 * w) + 273.15  # K

rho_Twp(T, w, p) = p / (287.05 * T * (1 + 1.6078 * w))  # kg_da/m3
rho_p(p) = p * rho_air_ref / p_air_ref

w_TRHp(T, RH, p) = 0.621945 * (RH * p_saturation(T)) / (p - (RH * p_saturation(T)))  # kg_w/kg_da

# Alduchov, O. A., and R. E. Eskridge, 1996: Improved Magnus Form Approximation of Saturation Vapor Pressure. J. Appl. Meteor. Climatol., 35, 601-609, https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2.
p_saturation(T) = 610.94 * exp(17.625 * (T - 273.15) / (T - 273.15 + 243.04))  # Pa


# Stull, R., 2011: Wet-Bulb Temperature from Relative Humidity and Air Temperature. J. Appl. Meteor. Climatol., 50, 2267-2269,  https://doi.org/10.1175/JAMC-D-11-0143.1.
wetbulb_TRH(T, RH) = 273.15 + (T - 273.15) * atan(0.151977 * sqrt(100 * RH + 8.313659)) +
    atan((T - 273.15) + 100 * RH) - atan(100 * RH - 1.676331) +
    0.00391838 * (100 * RH)^(3 / 2) * atan(0.023101 * 100 * RH) - 4.686035  # K

w_sat(T) = 0.621945 * p_saturation(T) / (patm - p_saturation(T))  # kg_w/kg_da

T_dew(w) = let pw = w * p_air_ref / (0.621945 + w)
    243.04 * log(pw / 610.94) / (17.625 - log(pw / 610.94)) + 273.15
end  # K
