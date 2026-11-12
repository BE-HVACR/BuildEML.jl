# Moist-air and liquid-water well-mixed control-volume nodes.
#
# Air (N-port, N >= 1, port flow direction not fixed):
#   AirMixingVolumeNodeN  — tracks T, w, p; set if_moisture_input=true to expose water_mflow (liquid water injection) and w_out ports
#
# Water (N-port, N >= 1, port flow direction not fixed):
#   WaterMixingVolumeNodeN — tracks T; incompressible, so no pressure state
#
# If added liquid water (water_mflow) also carries energy, inject that energy via heatport at the instantiation site.

using ModelingToolkitStandardLibrary.Thermal: HeatPort
using ModelingToolkitStandardLibrary.Blocks: Constant, RealInput, RealOutput


@mtkmodel AirMixingVolumeNodeN begin
    @structural_parameters begin
        N::Int = 3
        if_steady_state::Bool = false
        if_moisture_input::Bool = false
    end

    @parameters begin
        V = 1.0               # [m^3]
        Qflow_const = 0.0     # [W]
        δ = 1e-6
    end

    @components begin
        airport = [AirPort() for i in 1:N]

        heatport = HeatPort()
        self_heatport = HeatPort()
        pretem = PrescribedTemperature()
        heatflowsen = HeatFlowSensor()
        water_mflow_ = RealInput()
        if if_moisture_input
            water_mflow = RealInput()
            w_out = RealOutput()
        end
    end

    @variables begin
        T_mix(t), [guess = 293.15]
        w_mix(t), [guess = 0.008]
        p_mix(t), [guess = 101325.0]
    end

    @equations begin
        [airport[i].p ~ p_mix for i in 1:N]...

        [airport[i].T_ofo ~ T_mix for i in 1:N]...
        [airport[i].w_ofo ~ w_mix for i in 1:N]...

        pretem.T.u ~ T_mix
        connect(self_heatport, pretem.port)
        connect(heatport, heatflowsen.port_a)
        connect(self_heatport, heatflowsen.port_b)

        if !if_steady_state
            V * (rho_p(p_mix) / p_mix) * D(p_mix) ~ sum(airport[i].mflow for i in 1:N)
            (rho_p(p_mix) * V) * D(w_mix) ~
                sum(smooth_max(airport[i].mflow, 0.0; δ = δ) * (instream(airport[i].w_ofo) - w_mix) for i in 1:N) +
                water_mflow_.u
            (rho_p(p_mix) * V) * (cp_moistair_w(w_mix) * D(T_mix) + h_g(T_mix) * D(w_mix)) ~
                sum(smooth_max(airport[i].mflow, 0.0; δ = δ) *
                    (h_Tw(instream(airport[i].T_ofo), instream(airport[i].w_ofo)) - h_Tw(T_mix, w_mix)) for i in 1:N) +
                Qflow_const + heatflowsen.Q_flow.u
        else
            0 ~ sum(airport[i].mflow for i in 1:N)
            0 ~ sum(smooth_max(airport[i].mflow, 0.0; δ = δ) * (instream(airport[i].w_ofo) - w_mix) for i in 1:N) +
                water_mflow_.u
            0 ~ sum(smooth_max(airport[i].mflow, 0.0; δ = δ) *
                    (h_Tw(instream(airport[i].T_ofo), instream(airport[i].w_ofo)) - h_Tw(T_mix, w_mix)) for i in 1:N) +
                Qflow_const + heatflowsen.Q_flow.u
        end

        if if_moisture_input
            connect(water_mflow, water_mflow_)
            w_out.u ~ w_mix
        else
            water_mflow_.u ~ 0.0
        end
    end
end


@mtkmodel WaterMixingVolumeNodeN begin
    @structural_parameters begin
        N::Int = 2
        if_steady_state::Bool = false
    end

    @parameters begin
        V = 0.01          # [m^3]
    end

    @components begin
        ports = [WaterPort() for i in 1:N]

        heatport = HeatPort()
        self_heatport = HeatPort()
        pretem = PrescribedTemperature()
        heatflowsen = HeatFlowSensor()
    end

    @variables begin
        T_mix(t), [guess = 313.15]
        p_mix(t), [guess = 1.5e5]
    end

    @equations begin
        [ports[i].p ~ p_mix for i in 1:N]...
        0 ~ sum(ports[i].mflow for i in 1:N)

        [ports[i].T_ofo ~ T_mix for i in 1:N]...

        pretem.T.u ~ T_mix
        connect(self_heatport, pretem.port)
        connect(heatport, heatflowsen.port_a)
        connect(self_heatport, heatflowsen.port_b)

        if !if_steady_state
            (rho_water * V * cp_water) * D(T_mix) ~
                sum(smooth_max(ports[i].mflow, 0.0) * cp_water * (instream(ports[i].T_ofo) - T_mix) for i in 1:N) +
                heatflowsen.Q_flow.u
        else
            0 ~ sum(smooth_max(ports[i].mflow, 0.0) * cp_water * (instream(ports[i].T_ofo) - T_mix) for i in 1:N) +
                heatflowsen.Q_flow.u
        end
    end
end

