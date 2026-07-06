const _AIR_MEDIA = :air
const _WATER_MEDIA = :water

# MTK may emit bare `air`/`water` names for Symbol structural parameters.
const air = _AIR_MEDIA
const water = _WATER_MEDIA


"""
Ideal electric power model for heating mode using the core
Buildings.Fluid.HeatPumps.Carnot_TCon COP calculation.

Aligned with Carnot_TCon:
- COP = etaCarnot * COPCarnot * etaPL.
- TConAct = T_hot_ref + yPL*TAppCon_nominal.
- TEvaAct = T_cold_ref - yPL*TAppEva_nominal.
- yPL = max(Q_thermal, 0)/QCon_flow_nominal, limited to [yPL_min, yPL_max].
  Negative `Q_thermal` is ignored by this heating block.
- QCon_flow_nominal must be supplied by the user and must be positive, matching
  Carnot_TCon.
- etaPL is a polynomial in yPL: etaPL_a0 + etaPL_a1*yPL + ...
- If use_eta_Carnot_nominal=false, etaCarnot is computed from COP_nominal and
  the nominal temperatures, as in the Modelica model.

Differences from Carnot_TCon:
- This is a steady algebraic power block, not a four-port fluid model.
- It does not compute condenser or evaporator heat transfer from mass flow and
  leaving temperatures; those are supplied as `Q_thermal`, `T_hot_ref`
  (condenser side), and `T_cold_ref` (evaporator side).
- The part-load polynomial is fixed to cubic coefficients instead of a
  variable-length Modelica array.
- Modelica derives the default approach temperatures from the fluid media
  specific heat. This block uses `con_side_media` and `eva_side_media` instead:
  `:air` gives 5 K and `:water` gives 2 K.
- This is the heating-mode block. Connect `T_hot_ref` to the supply-air
  condenser-side temperature and `T_cold_ref` to the outdoor evaporator-side
  temperature.
  Use `HPCooling_Power_TEva` for cooling mode.

NOTE:
- Known limitation: At very small temperature lift, the Carnot COP can become
  very large, which makes the computed electric power very small. This is kept
  as an approximation for now; this block does not model compressor cycling,
  minimum pressure ratio, or other real-equipment limits.
"""
@mtkmodel HPHeating_Power_TCon begin
    @structural_parameters begin
        use_eta_Carnot_nominal::Bool = true
        con_side_media::Symbol = _AIR_MEDIA
        eva_side_media::Symbol = _AIR_MEDIA
    end

    @parameters begin
        # Carnot effectiveness and nominal fallback, matching Carnot_TCon.
        etaCarnot_nominal = 0.3
        TCon_nominal = 303.15
        TEva_nominal = 278.15

        # Approach temperatures, matching Carnot_TCon.
        # [K] condenser refrigerant temp above condenser-side reference temp
        TAppCon_nominal =
            con_side_media == _AIR_MEDIA ? 5.0 :
            con_side_media == _WATER_MEDIA ? 2.0 :
            error("Unsupported con_side_media $(con_side_media). Use :air or :water.")
        # [K] evaporator refrigerant temp below evaporator-side reference temp
        TAppEva_nominal =
            eva_side_media == _AIR_MEDIA ? 5.0 :
            eva_side_media == _WATER_MEDIA ? 2.0 :
            error("Unsupported eva_side_media $(eva_side_media). Use :air or :water.")

        COP_nominal = etaCarnot_nominal * (TCon_nominal + TAppCon_nominal) /
                      ((TCon_nominal + TAppCon_nominal) - (TEva_nominal - TAppEva_nominal))

        # Part-load normalization and limits. Must be positive, as in Carnot_TCon.
        QCon_flow_nominal
        yPL_min = 0.0
        yPL_max = 1.0

        # Part-load efficiency polynomial. The default etaPL is 1 at all loads.
        etaPL_a0 = 1.0
        etaPL_a1 = 0.0
        etaPL_a2 = 0.0
        etaPL_a3 = 0.0

        # Minimums for numerical protection.
        dT_min = 1.0
        etaPL_min = 1e-3

        # Minimum COP to avoid division by very small/negative values
        COP_heat_min = 1.0
        COP_min      = 0.5
    end

    @components begin
        Q_thermal  = RealInput()    # [W] heating load, positive in normal use
        T_hot_ref  = RealInput()    # [K] condenser-side leaving/reference temp
        T_cold_ref = RealInput()    # [K] evaporator-side leaving/reference temp

        power   = RealOutput()      # [W] electric power
        COP     = RealOutput()      # [-] effective COP
        dT      = RealOutput()      # [K] effective lift
        TConAct = RealOutput()      # [K] refrigerant-side condenser temperature
        TEvaAct = RealOutput()      # [K] refrigerant-side evaporator temperature
        yPL     = RealOutput()      # [-] part-load ratio
        etaPL   = RealOutput()      # [-] part-load efficiency multiplier
    end

    @variables begin
        Q_flow_(t)
        yPL_(t)
        etaPL_(t)
        etaCarnot_(t)
        TConAct_(t)
        TEvaAct_(t)
        dT_(t)
        COPCar_(t)
        COP_eff_(t)
        power_(t)
    end

    @equations begin
        Q_flow_ ~ smooth_max(Q_thermal.u, 0.0)
        yPL_ ~ smooth_clamp(Q_flow_ / QCon_flow_nominal, yPL_min, yPL_max)

        etaPL_ ~ smooth_max(
            etaPL_a0 + etaPL_a1*yPL_ + etaPL_a2*yPL_^2 + etaPL_a3*yPL_^3,
            etaPL_min
        )

        if use_eta_Carnot_nominal
            etaCarnot_ ~ etaCarnot_nominal
        else
            etaCarnot_ ~ COP_nominal /
                         ((TCon_nominal + TAppCon_nominal) /
                          ((TCon_nominal + TAppCon_nominal) - (TEva_nominal - TAppEva_nominal)))
        end

        TConAct_ ~ T_hot_ref.u + yPL_ * TAppCon_nominal
        TEvaAct_ ~ T_cold_ref.u - yPL_ * TAppEva_nominal
        dT_ ~ smooth_max(TConAct_ - TEvaAct_, dT_min)

        COPCar_ ~ TConAct_ / dT_
        COP_eff_ ~ smooth_max(smooth_max(etaCarnot_ * COPCar_ * etaPL_, COP_heat_min), COP_min)

        power_ ~ Q_flow_ / COP_eff_

        power.u   ~ power_
        COP.u     ~ COP_eff_
        dT.u      ~ dT_
        TConAct.u ~ TConAct_
        TEvaAct.u ~ TEvaAct_
        yPL.u     ~ yPL_
        etaPL.u   ~ etaPL_
    end
end


"""
Ideal electric power model for cooling mode using the core
Buildings.Fluid.Chillers.Carnot_TEva COP calculation.

Aligned with Carnot_TEva:
- COP = etaCarnot * COPCarnot * etaPL.
- COPCarnot uses TEvaAct as the useful temperature.
- TConAct = T_hot_ref + yPL*TAppCon_nominal.
- TEvaAct = T_cold_ref - yPL*TAppEva_nominal.
- yPL = max(-Q_thermal, 0)/(-QEva_flow_nominal), limited to
  [yPL_min, yPL_max]. Positive `Q_thermal` is ignored by this cooling block.
- QEva_flow_nominal must be supplied by the user and must be negative, matching
  Carnot_TEva.
- etaPL is a polynomial in yPL: etaPL_a0 + etaPL_a1*yPL + ...
- If use_eta_Carnot_nominal=false, etaCarnot is computed from COP_nominal and
  the nominal temperatures, as in the Modelica model.

Differences from Carnot_TEva:
- This is a steady algebraic power block, not a four-port fluid model.
- It does not compute condenser or evaporator heat transfer from mass flow and
  leaving temperatures; those are supplied as `Q_thermal`, `T_hot_ref`
  (condenser side), and `T_cold_ref` (evaporator side).
- The part-load polynomial is fixed to cubic coefficients instead of a
  variable-length Modelica array.
- Modelica derives the default approach temperatures from the fluid media
  specific heat. This block uses `con_side_media` and `eva_side_media` instead:
  `:air` gives 5 K and `:water` gives 2 K.
- This is the cooling-mode block. Connect `T_hot_ref` to the outdoor
  condenser-side temperature and `T_cold_ref` to the supply-air evaporator-side
  temperature. Use `HPHeating_Power_TCon` for heating mode.

NOTE:
- Known limitation: At very small temperature lift, the Carnot COP can become
  very large, which makes the computed electric power very small. This is kept
  as an approximation for now; this block does not model compressor cycling,
  minimum pressure ratio, or other real-equipment limits.
"""
@mtkmodel HPCooling_Power_TEva begin
    @structural_parameters begin
        use_eta_Carnot_nominal::Bool = true
        con_side_media::Symbol = _AIR_MEDIA
        eva_side_media::Symbol = _AIR_MEDIA
    end

    @parameters begin
        # Carnot effectiveness and nominal fallback, matching Carnot_TEva.
        etaCarnot_nominal = 0.3
        TCon_nominal = 303.15
        TEva_nominal = 278.15

        # Approach temperatures, matching Carnot_TEva.
        # [K] condenser refrigerant temp above condenser-side reference temp
        TAppCon_nominal =
            con_side_media == _AIR_MEDIA ? 5.0 :
            con_side_media == _WATER_MEDIA ? 2.0 :
            error("Unsupported con_side_media $(con_side_media). Use :air or :water.")
        # [K] evaporator refrigerant temp below evaporator-side reference temp
        TAppEva_nominal =
            eva_side_media == _AIR_MEDIA ? 5.0 :
            eva_side_media == _WATER_MEDIA ? 2.0 :
            error("Unsupported eva_side_media $(eva_side_media). Use :air or :water.")

        COP_nominal = etaCarnot_nominal * (TEva_nominal - TAppEva_nominal) /
                      ((TCon_nominal + TAppCon_nominal) - (TEva_nominal - TAppEva_nominal))

        # Part-load normalization and limits. Must be negative, as in Carnot_TEva.
        QEva_flow_nominal
        yPL_min = 0.0
        yPL_max = 1.0

        # Part-load efficiency polynomial. The default etaPL is 1 at all loads.
        etaPL_a0 = 1.0
        etaPL_a1 = 0.0
        etaPL_a2 = 0.0
        etaPL_a3 = 0.0

        # Minimums for numerical protection.
        dT_min = 1.0
        etaPL_min = 1e-3

        # Minimum COP to avoid division by very small/negative values
        COP_cool_min = 0.5
        COP_min      = 0.5
    end

    @components begin
        Q_thermal  = RealInput()    # [W] cooling load, negative in normal use
        T_hot_ref  = RealInput()    # [K] condenser-side leaving/reference temp
        T_cold_ref = RealInput()    # [K] evaporator-side leaving/reference temp

        power   = RealOutput()      # [W] electric power
        COP     = RealOutput()      # [-] effective cooling COP
        dT      = RealOutput()      # [K] effective lift
        TConAct = RealOutput()      # [K] refrigerant-side condenser temperature
        TEvaAct = RealOutput()      # [K] refrigerant-side evaporator temperature
        yPL     = RealOutput()      # [-] part-load ratio
        etaPL   = RealOutput()      # [-] part-load efficiency multiplier
    end

    @variables begin
        Q_flow_(t)
        yPL_(t)
        etaPL_(t)
        etaCarnot_(t)
        TConAct_(t)
        TEvaAct_(t)
        dT_(t)
        COPCar_(t)
        COP_eff_(t)
        power_(t)
    end

    @equations begin
        Q_flow_ ~ smooth_max(-Q_thermal.u, 0.0)
        yPL_ ~ smooth_clamp(Q_flow_ / (-QEva_flow_nominal), yPL_min, yPL_max)

        etaPL_ ~ smooth_max(
            etaPL_a0 + etaPL_a1*yPL_ + etaPL_a2*yPL_^2 + etaPL_a3*yPL_^3,
            etaPL_min
        )

        if use_eta_Carnot_nominal
            etaCarnot_ ~ etaCarnot_nominal
        else
            etaCarnot_ ~ COP_nominal /
                         ((TEva_nominal - TAppEva_nominal) /
                          ((TCon_nominal + TAppCon_nominal) - (TEva_nominal - TAppEva_nominal)))
        end

        TConAct_ ~ T_hot_ref.u + yPL_ * TAppCon_nominal
        TEvaAct_ ~ T_cold_ref.u - yPL_ * TAppEva_nominal
        dT_ ~ smooth_max(TConAct_ - TEvaAct_, dT_min)

        COPCar_ ~ TEvaAct_ / dT_
        COP_eff_ ~ smooth_max(smooth_max(etaCarnot_ * COPCar_ * etaPL_, COP_cool_min), COP_min)

        power_ ~ Q_flow_ / COP_eff_

        power.u   ~ power_
        COP.u     ~ COP_eff_
        dT.u      ~ dT_
        TConAct.u ~ TConAct_
        TEvaAct.u ~ TEvaAct_
        yPL.u     ~ yPL_
        etaPL.u   ~ etaPL_
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
