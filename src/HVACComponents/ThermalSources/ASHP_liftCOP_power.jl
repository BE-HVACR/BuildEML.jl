"""
Ideal electric power model for air-source heat pumps based on a lift-dependent COP.

COP correlation (Staffell et al., Energy & Environ. Sci. 2012, 5, 9291):

    COP_heat(ΔT) = 6.81 − 0.121·ΔT + 0.000630·ΔT²   (fitted for 15–60 K)

ΔT is clamped to [dT_min, dT_max] before applying the polynomial.
Cooling COP is derived as COP_cool ≈ max(COP_heat − 1, COP_cool_min).
A smooth blend selects between heating and cooling COP based on the sign of Q_thermal.
Power = |Q_thermal| / COP_eff.
"""
@mtkmodel ASHP_LiftCOP_Power begin
    @parameters begin
        # Staffell ASHP quadratic coefficients (heating COP)
        a0 = 6.81
        a1 = -0.121
        a2 = 0.000630

        # Computational ΔT range [K] (original fit: 15–60 K)
        dT_min = 15.0
        dT_max = 60.0

        # Offsets: map reference temps to effective refrigerant-side temps
        T_hot_offset  =  8.0   # [K] condenser hotter than outdoor air
        T_cold_offset = -2.0   # [K] evaporator colder than cold-side air

        # Minimum COPs to avoid division by very small/negative values
        COP_heat_min = 1.0
        COP_cool_min = 0.5
        COP_min      = 0.5
    end

    @components begin
        Q_thermal  = RealInput()    # [W] >0 heating, <0 cooling
        T_hot_ref  = RealInput()    # [K] hot-side reference temp
        T_cold_ref = RealInput()    # [K] cold-side reference temp

        power = RealOutput()        # [W] electric power
        COP   = RealOutput()        # [-] effective COP
        dT    = RealOutput()        # [K] effective lift
    end

    @variables begin
        dT_(t)
        COP_heat_(t)
        COP_cool_(t)
        COP_eff_(t)
        power_(t)
    end

    @equations begin
        dT_ ~ smooth_clamp((T_hot_ref.u + T_hot_offset) - (T_cold_ref.u + T_cold_offset), dT_min, dT_max)

        COP_heat_ ~ smooth_max(a0 + a1*dT_ + a2*dT_^2, COP_heat_min)

        # Cooling COP ≈ COP_heat - 1, with a floor
        COP_cool_ ~ smooth_max(COP_heat_ - 1.0, COP_cool_min)

        # Smooth blend: Q>0 → COP_heat, Q<0 → COP_cool
        COP_eff_ ~ smooth_max(soft_blend(COP_heat_, COP_cool_, Q_thermal.u), COP_min)

        power_ ~ smooth_abs(Q_thermal.u) / COP_eff_

        power.u ~ power_
        COP.u   ~ COP_eff_
        dT.u    ~ dT_
    end
end


"""Constant-COP electric power model: power = |Q_thermal| / COP_const."""
@mtkmodel ConstantCOP_Power begin
    @parameters begin
        COP_const = 1.0
    end

    @components begin
        Q_thermal = RealInput()    # [W] >0 heating, <0 cooling
        power     = RealOutput()   # [W] electric power
    end

    @equations begin
        power.u ~ smooth_abs(Q_thermal.u) / COP_const
    end
end
