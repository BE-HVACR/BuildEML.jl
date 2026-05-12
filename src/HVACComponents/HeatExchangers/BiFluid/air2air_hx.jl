"""
Constant-effectiveness sensible (and optional latent) air-to-air heat exchanger.

- Two air streams: side 1 (port_a1 → port_b1), side 2 (port_a2 → port_b2).
- Sensible: Q_sens = eff_sensible · Cmin · (T2_in − T1_in); positive means heat flows 2→1.
- Optional latent (if_latent=true): W_lat = eff_latent · Cwmin · (w2_in − w1_in).
- Pressure-drop modes: if_dp_quadratic takes precedence over if_dp_const; both default to false.
- Cp assumed constant; no flow reversal assumed (a→b is inlet→outlet).
"""
@mtkmodel Constant_HX begin
    @structural_parameters begin
        if_latent::Bool       = false
        if_dp_const::Bool     = false
        if_dp_quadratic::Bool = false
    end

    @parameters begin
        eff_sensible = 0.8        # [-]
        eff_latent   = 0.7        # [-]

        mflow_nominal = 0.2       # [kg/s]
        dp_nominal    = 50.0      # [Pa]
        Cflow_min     = 1e-4      # [W/K] regularization lower bound
    end

    @components begin
        port_a1 = AirPort()
        port_b1 = AirPort()
        port_a2 = AirPort()
        port_b2 = AirPort()
    end

    @variables begin
        m1_flow(t)
        m2_flow(t)

        T1_in(t)
        T1_out(t), [guess = 293.15]
        T2_in(t)
        T2_out(t), [guess = 293.15]

        w1_in(t)
        w1_out(t), [guess = 0.008]
        w2_in(t)
        w2_out(t), [guess = 0.008]

        C1_flow(t)
        C2_flow(t)
        Cmin_flow(t)
        C1_eff(t)
        C2_eff(t)

        Cw1_flow(t)
        Cw2_flow(t)
        Cwmin_flow(t)

        Q_sens(t)      # [W]
        W_lat(t)       # [kg/s]

        dp1(t)
        dp2(t)
    end

    @equations begin
        port_a1.mflow + port_b1.mflow ~ 0
        m1_flow ~ port_a1.mflow

        port_a2.mflow + port_b2.mflow ~ 0
        m2_flow ~ port_a2.mflow

        T1_in ~ instream(port_a1.T_ofo)
        T2_in ~ instream(port_a2.T_ofo)

        w1_in ~ instream(port_a1.w_ofo)
        w2_in ~ instream(port_a2.w_ofo)

        C1_flow   ~ smooth_abs(m1_flow) * cp_da
        C2_flow   ~ smooth_abs(m2_flow) * cp_da
        Cmin_flow ~ smooth_min(C1_flow, C2_flow)
        C1_eff    ~ smooth_max(C1_flow, Cflow_min)
        C2_eff    ~ smooth_max(C2_flow, Cflow_min)

        Q_sens ~ eff_sensible * Cmin_flow * (T2_in - T1_in)

        T1_out ~ T1_in + Q_sens / C1_eff
        T2_out ~ T2_in - Q_sens / C2_eff

        if !if_latent
            w1_out ~ w1_in
            w2_out ~ w2_in
            W_lat  ~ 0.0
        else
            Cw1_flow   ~ smooth_abs(m1_flow)
            Cw2_flow   ~ smooth_abs(m2_flow)
            Cwmin_flow ~ Cw1_flow + Cw2_flow - smooth_max(Cw1_flow, Cw2_flow)

            W_lat ~ eff_latent * Cwmin_flow * (w2_in - w1_in)

            Cw1_flow * (w1_out - w1_in) ~ W_lat
            Cw2_flow * (w2_out - w2_in) ~ -W_lat
        end

        port_b1.T_ofo ~ T1_out
        port_b2.T_ofo ~ T2_out

        port_b1.w_ofo ~ w1_out
        port_b2.w_ofo ~ w2_out

        if if_dp_quadratic
            dp1 ~ dp_nominal *
                   (m1_flow * smooth_abs(m1_flow)) /
                   (mflow_nominal * smooth_abs(mflow_nominal) + 1e-6)

            dp2 ~ dp_nominal *
                   (m2_flow * smooth_abs(m2_flow)) /
                   (mflow_nominal * smooth_abs(mflow_nominal) + 1e-6)

            port_b1.p ~ port_a1.p - dp1
            port_b2.p ~ port_a2.p - dp2
        elseif if_dp_const
            dp1 ~ dp_nominal
            dp2 ~ dp_nominal

            port_b1.p ~ port_a1.p - dp1
            port_b2.p ~ port_a2.p - dp2
        else
            port_b1.p ~ port_a1.p
            port_b2.p ~ port_a2.p
        end
    end
end
