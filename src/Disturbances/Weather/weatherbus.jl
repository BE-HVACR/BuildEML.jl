using ModelingToolkitStandardLibrary.Blocks
using DataInterpolations


function _right_periodic_pad(time_data::AbstractVector, y_data::AbstractVector, n_pad::Int)
    n_pad <= 0 && return time_data, y_data
    length(time_data) >= 2 || return time_data, y_data

    dt = time_data[end] - time_data[end - 1]
    dt > 0 || error("WeatherBus: time data must be strictly increasing for periodic padding.")

    # ReadEPW-aligned weather data has a t=0 anchor row, so use row 2 as
    # the first actual next-year sample when padding the right boundary.
    nper = length(time_data)
    pad_idx = [mod1(j + 1, nper) for j in 1:n_pad]
    time_pad = [time_data[end] + j * dt for j in 1:n_pad]

    return vcat(time_data, time_pad), vcat(y_data, y_data[pad_idx])
end


"""
    WeatherBus(df; name, time_col, interp_method, use_constant_pressure,
               radiation_time_shift_s, periodic_padding_steps)

Build a WeatherBus ODESystem from a processed weather DataFrame (e.g. from `ReadEPW`).
Exposes thermo, radiation, and wind signals as `RealOutput` ports.
Uses `AkimaInterpolation` by default; pass `interp_method` to override.
With Akima interpolation, `periodic_padding_steps=1` appends one next-year anchor
to reduce right-boundary artifacts. This approximates the Modelica Buildings
last-two-point boundary extrapolation more simply.

`df` must contain: `:time` [s], `:TDryBul` [K], `:TDewPoi` [K], `:relHum` [-],
`:TWetBul` [K], `:HumRat` [kg/kg_da], `:HGloHor/:HDifHor/:HDirNor/:HHorIR` [W/m²],
`:albedo` [-], `:winDir` [rad], `:winSpe` [m/s].
When `use_constant_pressure=true`, `pAtm` is fixed at 101325 Pa. Otherwise, `df`
must contain `:pAtm` [Pa].

Solar radiation outputs are time-shifted by `radiation_time_shift_s` (default 1800 s)
to match Modelica ReaderTMY3 convention.
"""
function WeatherBus(df::DataFrame; name::Symbol = :WeatherBus,
                                   time_col::Symbol = :time,
                                   interp_method = AkimaInterpolation,
                                   use_constant_pressure::Bool = true,
                                   radiation_time_shift_s::Real = 1800.0,
                                   periodic_padding_steps::Int = 1)
    periodic_padding_steps >= 0 || error("WeatherBus: periodic_padding_steps must be nonnegative.")

    # time axis
    time_data    = Float64.(df[!, time_col])
    shift_s = Float64(radiation_time_shift_s)

    # thermo
    TDryBul_data = Float64.(df[!, :TDryBul])
    TDewPoi_data = Float64.(df[!, :TDewPoi])
    relHum_data  = Float64.(df[!, :relHum])
    pAtm_data    = use_constant_pressure ? fill(101325.0, length(time_data)) : Float64.(df[!, :pAtm])
    TWetBul_data = Float64.(df[!, :TWetBul])
    HumRat_data  = Float64.(df[!, :HumRat])

    # radiation
    HGloHor_data = Float64.(df[!, :HGloHor])
    HDifHor_data = Float64.(df[!, :HDifHor])
    HDirNor_data = Float64.(df[!, :HDirNor])
    HHorIR_data  = Float64.(df[!, :HHorIR])
    albedo_data  = Float64.(df[!, :albedo])

    # wind
    winDir_data  = Float64.(df[!, :winDir])
    winSpe_data  = Float64.(df[!, :winSpe])

    time_data_pad, TDryBul_data = _right_periodic_pad(time_data, TDryBul_data, periodic_padding_steps)
    _, TDewPoi_data = _right_periodic_pad(time_data, TDewPoi_data, periodic_padding_steps)
    _, relHum_data  = _right_periodic_pad(time_data, relHum_data,  periodic_padding_steps)
    _, pAtm_data    = _right_periodic_pad(time_data, pAtm_data,    periodic_padding_steps)
    _, TWetBul_data = _right_periodic_pad(time_data, TWetBul_data, periodic_padding_steps)
    _, HumRat_data  = _right_periodic_pad(time_data, HumRat_data,  periodic_padding_steps)
    _, HHorIR_data  = _right_periodic_pad(time_data, HHorIR_data,  periodic_padding_steps)
    _, albedo_data  = _right_periodic_pad(time_data, albedo_data,  periodic_padding_steps)
    _, winDir_data  = _right_periodic_pad(time_data, winDir_data,  periodic_padding_steps)
    _, winSpe_data  = _right_periodic_pad(time_data, winSpe_data,  periodic_padding_steps)

    # Add enough right padding to cover the solar-radiation time shift and
    # retain the Akima endpoint stabilization used by the other weather signals.
    rad_padding_steps = periodic_padding_steps
    if shift_s > 0 && length(time_data) >= 2
        dt = time_data[2] - time_data[1]
        rad_padding_steps += max(1, ceil(Int, shift_s / dt))
    end

    time_data_rad, HGloHor_data_rad = _right_periodic_pad(time_data, HGloHor_data, rad_padding_steps)
    _, HDifHor_data_rad = _right_periodic_pad(time_data, HDifHor_data, rad_padding_steps)
    _, HDirNor_data_rad = _right_periodic_pad(time_data, HDirNor_data, rad_padding_steps)
    time_data_rad = time_data_rad .- shift_s

    @named clk          = ContinuousClock()

    @named itp_TDryBul  = Interpolation(interp_method, TDryBul_data, time_data_pad)
    @named itp_TDewPoi  = Interpolation(interp_method, TDewPoi_data, time_data_pad)
    @named itp_relHum   = Interpolation(interp_method, relHum_data,  time_data_pad)
    @named itp_pAtm     = Interpolation(interp_method, pAtm_data,    time_data_pad)
    @named itp_TWetBul  = Interpolation(interp_method, TWetBul_data, time_data_pad)
    @named itp_HumRat   = Interpolation(interp_method, HumRat_data,  time_data_pad)
    @named itp_HGloHor  = Interpolation(interp_method, HGloHor_data_rad, time_data_rad)
    @named itp_HDifHor  = Interpolation(interp_method, HDifHor_data_rad, time_data_rad)
    @named itp_HDirNor  = Interpolation(interp_method, HDirNor_data_rad, time_data_rad)
    @named itp_HHorIR   = Interpolation(interp_method, HHorIR_data,  time_data_pad)
    @named itp_albedo   = Interpolation(interp_method, albedo_data,  time_data_pad)
    @named itp_winDir   = Interpolation(interp_method, winDir_data,  time_data_pad)
    @named itp_winSpe   = Interpolation(interp_method, winSpe_data,  time_data_pad)

    @named TDryBul = RealOutput()
    @named TDewPoi = RealOutput()
    @named relHum  = RealOutput()
    @named pAtm    = RealOutput()
    @named TWetBul = RealOutput()
    @named HumRat  = RealOutput()

    @named HGloHor = RealOutput()
    @named HDifHor = RealOutput()
    @named HDirNor = RealOutput()
    @named HHorIR  = RealOutput()
    @named albedo  = RealOutput()

    @named winDir  = RealOutput()
    @named winSpe  = RealOutput()

    eqs = [
        connect(clk.output, itp_TDryBul.input)
        connect(clk.output, itp_TDewPoi.input)
        connect(clk.output, itp_relHum.input)
        connect(clk.output, itp_pAtm.input)
        connect(clk.output, itp_TWetBul.input)
        connect(clk.output, itp_HumRat.input)

        connect(clk.output, itp_HGloHor.input)
        connect(clk.output, itp_HDifHor.input)
        connect(clk.output, itp_HDirNor.input)
        connect(clk.output, itp_HHorIR.input)
        connect(clk.output, itp_albedo.input)

        connect(clk.output, itp_winDir.input)
        connect(clk.output, itp_winSpe.input)

        connect(itp_TDryBul.output, TDryBul)
        connect(itp_TDewPoi.output, TDewPoi)
        connect(itp_relHum.output,  relHum)
        connect(itp_pAtm.output,    pAtm)
        connect(itp_TWetBul.output, TWetBul)
        connect(itp_HumRat.output,  HumRat)

        connect(itp_HGloHor.output, HGloHor)
        connect(itp_HDifHor.output, HDifHor)
        connect(itp_HDirNor.output, HDirNor)
        connect(itp_HHorIR.output,  HHorIR)
        connect(itp_albedo.output,  albedo)

        connect(itp_winDir.output,  winDir)
        connect(itp_winSpe.output,  winSpe)
    ]

    ODESystem(eqs, t; name,
        systems = [clk,
                   itp_TDryBul, itp_TDewPoi, itp_relHum, itp_pAtm, itp_TWetBul, itp_HumRat,
                   itp_HGloHor, itp_HDifHor, itp_HDirNor, itp_HHorIR, itp_albedo,
                   itp_winDir, itp_winSpe,
                   TDryBul, TDewPoi, relHum, pAtm, TWetBul, HumRat,
                   HGloHor, HDifHor, HDirNor, HHorIR, albedo,
                   winDir, winSpe])
end



"""
    WeatherBus(; name, df, time_col, interp_method, use_constant_pressure,
               radiation_time_shift_s, periodic_padding_steps, kwargs...)

Keyword-friendly wrapper for MTK compatibility: `@named wea = WeatherBus(df_weather)`
works after MTK rewrites the call to `WeatherBus(; name=:wea, df_weather=df_weather)`.
Picks the DataFrame from `df` keyword or, if absent, the first DataFrame-valued kwarg.
"""
function WeatherBus(; name::Symbol = :WeatherBus,
                      time_col::Symbol = :time,
                      df::Union{Nothing,AbstractDataFrame} = nothing,
                      interp_method = AkimaInterpolation,
                      use_constant_pressure::Bool = true,
                      radiation_time_shift_s::Real = 1800.0,
                      periodic_padding_steps::Int = 1,
                      kwargs...)
    df_ = df

    if df_ === nothing
        for (_, v) in kwargs
            if v isa AbstractDataFrame
                df_ = v
                break
            end
        end
    end

    if df_ === nothing
        keylist = isempty(kwargs) ? "(none)" : join(string.(keys(kwargs)), ", ")
        error("WeatherBus: missing DataFrame. "*
              "Pass a positional `WeatherBus(df)` or keyword `WeatherBus(df=...)`. "*
              "Got kwargs: $keylist")
    end

    return WeatherBus(df_; name=name, time_col=time_col, interp_method=interp_method,
                      use_constant_pressure=use_constant_pressure,
                      radiation_time_shift_s=radiation_time_shift_s,
                      periodic_padding_steps=periodic_padding_steps)
end
