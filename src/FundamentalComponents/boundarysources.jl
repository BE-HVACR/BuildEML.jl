"""
Unified wrapper for common 1-port air boundary conditions.

- `role = :source` with `mode = :mflow` dispatches to `AirSource_mflow`
- `role = :source` with `mode = :pressure` dispatches to `AirSource_pT`
- `role = :sink`   with `mode = :mflow` dispatches to `AirSink_mflow`
- `role = :sink`   with `mode = :pressure` dispatches to `AirSink_p`

This wrapper unifies the API entry point
"""
function AirBoundary(;
    name::Symbol = :AirBoundary,
    role::Symbol = :source,
    mode::Symbol = :pressure,
    use_w_in::Bool = false,
    p_par::Real = 101325.0,
    w_par::Real = 0.008,
    T_backflow::Real = 293.15,
    w_backflow::Real = 0.008,
)
    if role == :source
        if mode == :mflow
            return AirSource_mflow(; name, use_w_in, w_par)
        elseif mode == :pressure
            return AirSource_pT(; name, use_w_in, p_par, w_par)
        else
            throw(ArgumentError("Unsupported AirBoundary mode=$(mode) for role=:source. Use :mflow or :pressure."))
        end
    elseif role == :sink
        if mode == :mflow
            return AirSink_mflow(; name, T_backflow, w_backflow)
        elseif mode == :pressure
            return AirSink_p(; name, T_backflow, w_backflow)
        else
            throw(ArgumentError("Unsupported AirBoundary mode=$(mode) for role=:sink. Use :mflow or :pressure."))
        end
    else
        throw(ArgumentError("Unsupported AirBoundary role=$(role). Use :source or :sink."))
    end
end


"""
mflow.u should be input as positive
it will be automatically converted into negative in this function
"""
@mtkmodel AirSource_mflow begin
    @structural_parameters begin
        use_w_in::Bool=false
    end

    @parameters begin
        w_par::Real=0.008
    end

    @components begin
        port = AirPort()

        mflow = RealInput()
        T_in = RealInput()
        w_in = RealInput()
    end

    @equations begin
        port.mflow ~ -mflow.u
        
        port.T_ofo ~ T_in.u

        if use_w_in
            port.w_ofo ~ w_in.u
        else
            port.w_ofo ~ w_par
        end
    end
end



@mtkmodel AirSource_pT begin
    @structural_parameters begin
        use_w_in::Bool=false
    end

    @parameters begin
        p_par::Real=101325.0
        w_par::Real=0.008
    end

    @components begin
        port = AirPort()

        T_in = RealInput()
        w_in = RealInput()
    end

    @equations begin
        port.p ~ p_par
        port.T_ofo ~ T_in.u

        if use_w_in
            port.w_ofo ~ w_in.u
        else
            port.w_ofo ~ w_par
        end
    end
end


@mtkmodel AirSink_mflow begin
    @parameters begin
        T_backflow = 293.15 
        w_backflow = 0.008
    end

    @components begin
        port = AirPort()

        mflow = RealInput()

        T_exh = RealOutput()
        w_exh = RealOutput()
    end

    @equations begin
        T_exh.u ~ instream(port.T_ofo)
        w_exh.u ~ instream(port.w_ofo)
        port.mflow ~ mflow.u

        port.T_ofo ~ T_backflow
        port.w_ofo ~ w_backflow
    end
end


@mtkmodel AirSink_p begin
    @parameters begin
        T_backflow = 293.15 
        w_backflow = 0.008
    end

    @components begin
        port = AirPort()

        p = RealInput()

        T_exh = RealOutput()
        w_exh = RealOutput()
        mflow_exh = RealOutput()
    end

    @equations begin
        T_exh.u ~ instream(port.T_ofo)
        w_exh.u ~ instream(port.w_ofo)
        mflow_exh.u ~ port.mflow
        port.p ~ p.u

        port.T_ofo ~ T_backflow
        port.w_ofo ~ w_backflow
    end
end




@mtkmodel WaterSource_mflow begin
    @components begin
        port = WaterPort()

        mflow = RealInput()
        T_in = RealInput()
    end

    @equations begin
        port.mflow ~ -mflow.u
      
        port.T_ofo ~ T_in.u
    end
end


@mtkmodel WaterSink_pT begin
    @parameters begin
        p_back = 1e5
        T_back = 293.15
    end
    @components begin
        port = WaterPort()
    end
    @equations begin
        port.p     ~ p_back
        port.T_ofo ~ T_back
    end
end


@mtkmodel WaterBoundaryNode2_p begin
    @structural_parameters begin
        use_p_in::Bool   = false  
    end

    @parameters begin
        p_default = 1.5e5          # [Pa] 
    end

    @components begin
        port_a = WaterPort()
        port_b = WaterPort()

        p_in = RealInput()
    end

    @variables begin
        p_boundary(t)                              
    end

    @equations begin
        if use_p_in
            p_boundary ~ p_in.u
        else
            p_boundary ~ p_default
        end
        
        port_a.p ~ p_boundary
        port_b.p ~ p_boundary
        port_b.T_ofo ~ instream(port_a.T_ofo) 
        port_a.mflow + port_b.mflow ~ 0

    end
end

