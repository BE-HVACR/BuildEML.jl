"""
Variable-opening air damper with selectable opening characteristic.

Ports: `port_a` (inlet), `port_b` (outlet), `y` (opening signal [0, 1]).
Air state (T, w) passes through unchanged; only pressure is modified.

Two modes via `use_exponential`:
- false (default): power law, `f_y = y_leak + (1 - y_leak) * y^n_open`, the opening
  characteristic of MBL's two-way valves (`TwoWayLinear` at `n_open = 1`,
  `TwoWayQuickOpening` at `n_open = 1/alpha`) applied to an air damper.
- true: `Buildings.Fluid.Actuators.Dampers.Exponential` characteristic.

`dpFixed_nominal` lumps the fixed resistance of the rest of the branch in series
with the damper (both modes); set it to 0 for the damper alone.

Both modes regularize the flow below `deltaM * mflow_nominal` with the polynomial
transition of MBL's `basicFlowFunction_dp`; `deltaM` is MBL's `deltaM`, with the
same default (0.3).
"""
@mtkmodel Damper begin
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

        deltaM           = 0.3        # [-] turbulent transition at deltaM * mflow_nominal
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
        k_tot(t), [irreducible = true, guess = 0.01]  # [(kg·m)^0.5 Pa^-0.5]
        kThetaSqRt(t)  # [-]   sqrt of loss coefficient (exponential mode)
        dp_turb(t)     # [Pa]  turbulent transition pressure difference
        dp_norm(t)     # [-]   dp / dp_turb
    end

    @equations begin
        port_b.T_ofo ~ instream(port_a.T_ofo)
        port_b.w_ofo ~ instream(port_a.w_ofo)
        port_a.T_ofo ~ port_a.T_guess
        port_a.w_ofo ~ port_a.w_guess

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
            k_tot ~ mflow_nominal * f_y /
                    sqrt(smooth_max(dpDamper_nominal + dpFixed_nominal, dp_eps))
            kThetaSqRt ~ 0.0
        end

        dp ~ port_a.p - port_b.p
        sgn_dp ~ dp / (smooth_abs(dp) + dp_eps)

        dp_turb ~ (deltaM * mflow_nominal / smooth_max(k_tot, 1e-12))^2
        dp_norm ~ dp / dp_turb
        m_flow ~ ifelse(
            smooth_abs(dp) * smooth_max(k_tot, 1e-12)^2 > (deltaM * mflow_nominal)^2,
            k_tot * sqrt(smooth_abs(dp)) * sgn_dp,
            (1.40625 + (0.15625 * dp_norm^2 - 0.5625) * dp_norm^2) *
                (deltaM * mflow_nominal) * dp_norm,
        )

        port_a.mflow ~ m_flow
        port_a.mflow + port_b.mflow ~ 0
    end
end
