"""
Variable-opening air damper with pressure-drop/flow characteristic.

Ports: `port_a` (inlet), `port_b` (outlet), `y` (opening signal [0, 1]).
Air state (T, w) passes through unchanged; only pressure is modified.

Two modes via `use_exponential`:
- false (default): power-law, `f_y = y_leak + (1 - y_leak) * y^n_open`
- true: MBL exponentialDamper characteristic (piecewise polynomial/exponential,
  breakpoints yL/yU); optional series resistance `dpFixed_nominal` combined in series.
"""
@mtkmodel DamperExponential begin
    @structural_parameters begin
        from_dp::Bool                = true   # only from_dp = true supported
        linearized::Bool             = false  # placeholder
        homotopyInitialization::Bool = false  # placeholder
        use_exponential::Bool        = false
    end

    @parameters begin
        mflow_nominal    = 0.1        # [kg/s]
        dpDamper_nominal = 10.0       # [Pa]
        dpFixed_nominal  = 190.0      # [Pa]

        rho_nominal      = 1.2        # [kg/m^3] reserved, not used yet

        y_leak           = 0.05       # [-] leakage factor at y=0
        n_open           = 2.0        # [-] power-law opening exponent
        a                = -1.51      # [-] exponential damper coefficient
        b                = 9.45       # [-] exponential damper coefficient
        yL               = 15 / 90   # [-] lower breakpoint
        yU               = 55 / 90   # [-] upper breakpoint
        k1               = 0.45      # [-] loss coefficient at y=1
        l                = 1e-4      # [-] leakage ratio k(y=0)/k(y=1)

        dp_eps           = 1e-3       # [Pa] regularization
    end

    @components begin
        port_a = AirPort()
        port_b = AirPort()
        y      = RealInput()          # opening signal [0, 1]
    end

    @variables begin
        m_flow(t)      # [kg/s]
        dp(t)          # [Pa]
        f_y(t)         # [-]   effective opening factor
        y_eff(t)       # [-]   clipped opening command
        sgn_dp(t)      # [-]   sign of dp
        k_dam(t)       # [(kg·m)^0.5 Pa^-0.5]
        k_tot(t)       # [(kg·m)^0.5 Pa^-0.5]
        kThetaSqRt(t)  # [-]   sqrt of loss coefficient (exponential mode)
    end

    @equations begin
        port_b.T_ofo ~ instream(port_a.T_ofo)
        port_b.w_ofo ~ instream(port_a.w_ofo)

        y_eff ~ smooth_clamp(y.u, 0.0, 1.0)
        if use_exponential
            kThetaSqRt ~ ifelse(
                y_eff < yL,
                sqrt(exp(
                    ((log(k1 / smooth_max(l, 1e-12)^2) - b - a) / yL^2) * y_eff^2 +
                    ((-b * yL - 2 * log(k1 / smooth_max(l, 1e-12)^2) + 2 * b + 2 * a) / yL) * y_eff +
                    log(k1 / smooth_max(l, 1e-12)^2)
                )),
                ifelse(
                    y_eff > yU,
                    sqrt(exp(
                        ((log(k1) - a) / (yU^2 - 2 * yU + 1)) * y_eff^2 +
                        ((-b * yU^2 - 2 * log(k1) * yU - (-2 * b - 2 * a) * yU - b) / (yU^2 - 2 * yU + 1)) * y_eff +
                        ((log(k1) * yU^2 + b * yU^2 + (-2 * b - 2 * a) * yU + b + a) / (yU^2 - 2 * yU + 1))
                    )),
                    sqrt(exp(a + b * (1.0 - y_eff)))
                ),
            )
            # k_dam = m_flow/sqrt(dpDamper), consistent with PartialDamperExponential
            k_dam ~ (mflow_nominal / sqrt(smooth_max(dpDamper_nominal, dp_eps))) *
                    (sqrt(k1) / smooth_max(kThetaSqRt, 1e-12))
            k_tot ~ ifelse(
                dpFixed_nominal > dp_eps,
                sqrt(1.0 / (
                    1.0 / (mflow_nominal / sqrt(smooth_max(dpFixed_nominal, dp_eps)))^2 +
                    1.0 / smooth_max(k_dam^2, 1e-12)
                )),
                k_dam
            )
            f_y ~ 1.0
        else
            f_y ~ y_leak + (1.0 - y_leak) * y_eff^n_open
            k_dam ~ 0.0
            k_tot ~ 0.0
            kThetaSqRt ~ 0.0
        end

        dp ~ port_a.p - port_b.p
        sgn_dp ~ dp / (smooth_abs(dp) + dp_eps)

        if use_exponential
            m_flow ~ k_tot * sqrt(smooth_abs(dp)) * sgn_dp
        else
            m_flow ~ mflow_nominal *
                     f_y *
                     sqrt(smooth_abs(dp) / smooth_max(dpDamper_nominal + dpFixed_nominal, dp_eps)) *
                     sgn_dp
        end

        port_a.mflow ~ m_flow
        port_a.mflow + port_b.mflow ~ 0
    end
end
