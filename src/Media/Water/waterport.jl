"""
    WaterPort

Water-side connector carrying pressure, stream temperature, and mass flow rate.
"""
@connector WaterPort begin
    @parameters begin
        p_guess = 1.5e5        # [Pa]
        T_guess = 313.15       # [K]
        mflow_guess = 0.02     # [kg/s], mass flow rate
    end
    @variables begin
        p(t), [guess = p_guess]
        T_ofo(t), [guess = T_guess, connect = Stream]
        mflow(t), [guess = mflow_guess, connect = Flow]  # >0 = into component
    end
end
