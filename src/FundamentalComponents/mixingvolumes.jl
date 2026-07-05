# Moist-air and liquid-water well-mixed control-volume nodes.
#
# Air (N-port, arbitrary topology via port_role):
#   AirMixingVolumeNodeN  — tracks T, w, p; set if_moisture_input=true to expose water_mflow (liquid water injection) and w_out ports
#
# Water (2-port only):
#   WaterMixingVolumeNode2 — tracks T; port1 = inlet, port2 = outlet
#
# If added liquid water (water_mflow) also carries energy, inject that energy via heatport at the instantiation site.

using ModelingToolkitStandardLibrary.Thermal: HeatPort
using ModelingToolkitStandardLibrary.Blocks: Constant, RealInput, RealOutput


function _build_airmixingvolume_n_system(;
    name::Symbol,
    N::Int,
    if_steady_state::Bool = false,
    V::Real = 1.0,
    Qflow_const::Real = 0.0,
    delta::Real = 1e-6,
    if_moisture_input::Bool = false,
    port_role::Vector{Symbol} = [i == 1 ? :inlet : :outlet for i in 1:N],
)
    N >= 1 || throw(ArgumentError("AirMixingVolumeNodeN requires N >= 1."))
    length(port_role) == N || throw(ArgumentError("port_role must have length N."))
    all(role -> role in (:inlet, :outlet), port_role) ||
        throw(ArgumentError("port_role entries must be :inlet or :outlet."))

    ports = [AirPort(name = Symbol(:airport_, i)) for i in 1:N]

    @parameters begin
        Vp = V
        Qconst = Qflow_const
        δ = delta
    end

    @variables begin
        T_mix(t), [guess = 293.15]
        w_mix(t), [guess = 0.008]
        p_mix(t), [guess = 101325.0]
        port_T(t)[1:N], [guess = 293.15]
        port_w(t)[1:N], [guess = 0.008]
    end

    @named heatport = HeatPort()
    @named self_heatport = HeatPort()
    @named pretem = PrescribedTemperature()
    @named heatflowsen = HeatFlowSensor()
    @named water_mflow_ = RealInput()

    subs = ModelingToolkit.AbstractSystem[ports..., heatport, self_heatport, pretem, heatflowsen, water_mflow_]
    eqs = Equation[]

    append!(eqs, [ports[i].p ~ p_mix for i in 1:N])
    for i in 1:N
        if port_role[i] == :inlet
            push!(eqs, port_T[i] ~ instream(ports[i].T_ofo))
            push!(eqs, port_w[i] ~ instream(ports[i].w_ofo))
        else
            push!(eqs, port_T[i] ~ ports[i].T_ofo)
            push!(eqs, port_w[i] ~ ports[i].w_ofo)
        end
    end
    append!(eqs, [ports[i].T_ofo ~ T_mix for i in 1:N])
    append!(eqs, [ports[i].w_ofo ~ w_mix for i in 1:N])

    push!(eqs, pretem.T.u ~ T_mix)
    push!(eqs, connect(self_heatport, pretem.port))
    push!(eqs, connect(heatport, heatflowsen.port_a))
    push!(eqs, connect(self_heatport, heatflowsen.port_b))

    total_mflow = sum(port.mflow for port in ports)
    moisture_rhs = sum(ports[i].mflow * (port_w[i] - w_mix) for i in 1:N) + water_mflow_.u * (1.0 - w_mix)
    energy_rhs = sum(ports[i].mflow * (h_Tw(port_T[i], port_w[i]) - h_Tw(T_mix, w_mix)) for i in 1:N) +
                 Qconst + heatflowsen.Q_flow.u

    if !if_steady_state
        push!(eqs,
            Vp * ((-rho_p(p_mix) / T_mix) * D(T_mix) +
                  (-rho_p(p_mix) * 1.6078 / (1 + 1.6078 * w_mix)) * D(w_mix) +
                  (rho_p(p_mix) / p_mix) * D(p_mix)) ~ total_mflow
        )
        push!(eqs, (rho_p(p_mix) * Vp) * D(w_mix) ~ moisture_rhs)
        push!(eqs, (rho_p(p_mix) * Vp) * cp_moistair_w(w_mix) * D(T_mix) ~ energy_rhs)
    else
        push!(eqs, 0 ~ total_mflow)
        push!(eqs, 0 ~ moisture_rhs)
        push!(eqs, 0 ~ energy_rhs)
    end

    defaults = Dict()

    if if_moisture_input
        @named water_mflow = RealInput()
        @named w_out = RealOutput()
        push!(subs, water_mflow)
        push!(subs, w_out)
        push!(eqs, connect(water_mflow, water_mflow_))
        push!(eqs, w_out.u ~ w_mix)
    else
        push!(eqs, water_mflow_.u ~ 0.0)
    end

    return ODESystem(eqs, t; name, systems = subs, defaults = defaults)
end


function AirMixingVolumeNodeN(;
    name::Symbol = :AirMixingVolumeNodeN,
    N::Int = 3,
    if_steady_state::Bool = false,
    if_moisture_input::Bool = false,
    V::Real = 1.0,
    Qflow_const::Real = 0.0,
    delta::Real = 1e-6,
    port_role::Vector{Symbol} = [i == 1 ? :inlet : :outlet for i in 1:N],
)
    return _build_airmixingvolume_n_system(;
        name,
        N,
        if_steady_state,
        V,
        Qflow_const,
        delta,
        if_moisture_input,
        port_role,
    )
end









"""
Prefixed direction: port1 must be inlet and port2 must be outlet for this mixingvolume. 
"""
@mtkmodel WaterMixingVolumeNode2 begin
    @structural_parameters begin
        if_steady_state::Bool = false
    end

    @parameters begin
        V = 0.01          # [m^3]
    end

    @components begin
        port1 = WaterPort()
        port2 = WaterPort()

        heatport = HeatPort()
        self_heatport = HeatPort()
        pretem = PrescribedTemperature()
        heatflowsen = HeatFlowSensor()
    end

    @variables begin
        T_mix(t), [guess = 313.15]

        port1_T(t), [guess = 313.15]
        port2_T(t), [guess = 313.15]
    end

    @equations begin
        port1.p ~ port2.p
        0 ~ port1.mflow + port2.mflow

        port1.T_ofo ~ T_mix
        port2.T_ofo ~ T_mix

        port1_T ~ instream(port1.T_ofo)
        port2_T ~ port2.T_ofo

        pretem.T.u ~ T_mix
        connect(self_heatport, pretem.port)
        connect(heatport, heatflowsen.port_a)
        connect(self_heatport, heatflowsen.port_b)

        if !if_steady_state
            (rho_water * V * cp_water) * D(T_mix) ~
                port1.mflow * cp_water * port1_T +
                port2.mflow * cp_water * port2_T +
                heatflowsen.Q_flow.u
        else
            0 ~ port1.mflow * cp_water * port1_T +
                port2.mflow * cp_water * port2_T +
                heatflowsen.Q_flow.u
        end
    end
end
