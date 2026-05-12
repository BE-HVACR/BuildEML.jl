"""
    AirPort

Moist-air connector carrying pressure, outflow-only temperature, outflow-only humidity ratio, and dry air mass flow rate.
"""
@connector AirPort begin
    @parameters begin
        p_guess = 101325.0      # [Pa]
        T_guess = 293.15        # [K]
        w_guess = 0.008         # [kg_w/kg_da]
        mflow_guess = 0.1       # [kg/s], dry air mass flow rate
    end
    @variables begin
        p(t), [guess = p_guess]
        T_ofo(t), [guess = T_guess, connect = Stream]
        w_ofo(t), [guess = w_guess, connect = Stream]
        mflow(t), [guess = mflow_guess, connect = Flow]  # >0 = into component
    end
end
