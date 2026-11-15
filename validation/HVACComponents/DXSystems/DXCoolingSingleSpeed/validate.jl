# Validation of DXCoolingSingleSpeed against MBL v11
# Buildings.Fluid.DXSystems.Cooling.AirSource.Validation.SingleSpeed.
#
# SingleSpeed.mo is the MBL validation model; its input tables (coil inlet state,
# condenser inlet temperature, on/off) are parsed from it, and its source pressure
# pulse is reproduced. Buildings_Fluid_DXSystems_Cooling_AirSource_Validation_SingleSpeed.txt
# is the MBL Dymola reference result file (Resources/ReferenceResults/Dymola), unchanged.
#
# Differences from the MBL model:
# - accumulated-water reevaporation is not modeled (MBL computes it, but the air
#   flow is zero whenever the coil is off in this case)
# - input steps are replaced by 1 ms linear ramps

using BuildEML
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using DifferentialEquations
using DataInterpolations
using Plots
using Printf
using Statistics

include(joinpath(pkgdir(BuildEML), "validation", "helpers.jl"))

const DIR = @__DIR__

src = read(joinpath(DIR, "SingleSpeed.mo"), String)
const ITP_TC, kn1 = ramped(read_timetable(src, "TCIn")[1], read_timetable(src, "TCIn")[2] .+ 273.15)
const ITP_TE, kn2 = ramped(read_timetable(src, "TEvaIn")[1], read_timetable(src, "TEvaIn")[2] .+ 273.15)
const ITP_W, kn3 = ramped(read_timetable(src, "XEvaIn")...)   # humidity ratio, dry-air basis
const ITP_ON, kn4 = ramped(read_timetable(src, "plr_onOff")...)
# Pulse(offset = 101325, amplitude = 1141, startTime = 25200, period = 36000, width = 100 %, nperiod = 1)
const ITP_P, kn5 = ramped([0.0, 25200.0, 25200.0, 61200.0, 61200.0, 86400.0],
    [101325.0, 101325.0, 102466.0, 102466.0, 101325.0, 101325.0])

in_TC(x) = ITP_TC(x)
in_TE(x) = ITP_TE(x)
in_W(x) = ITP_W(x)
in_ON(x) = ITP_ON(x)
in_P(x) = ITP_P(x)
@register_symbolic in_TC(x)
@register_symbolic in_TE(x)
@register_symbolic in_W(x)
@register_symbolic in_ON(x)
@register_symbolic in_P(x)

@mtkmodel ValidationSource begin
    @components begin
        port = AirPort()
    end
    @equations begin
        port.p ~ in_P(t)
        port.T_ofo ~ in_TE(t)
        port.w_ofo ~ in_W(t)
    end
end

@mtkmodel SingleSpeedValidation begin
    @components begin
        sou = ValidationSource()
        sinSpeDX = DXCoolingSingleSpeed(
            perf_curves = DX_CURVE_II,
            Qflow_nominal = -10500.0,
            COP_nominal = 3.0,
            mflow_nominal = 1.72,
            SHR_nominal = 0.798655,
            dp_nominal = 1141.0,
        )
        sin = AirBoundary(role = :sink, mode = :pressure)
    end
    @equations begin
        connect(sou.port, sinSpeDX.port_a)
        connect(sinSpeDX.port_b, sin.port)
        sin.p.u ~ 101325.0
        sinSpeDX.on.u ~ in_ON(t)
        sinSpeDX.TConIn.u ~ in_TC(t)
    end
end

@mtkcompile sys = SingleSpeedValidation()

# MBL: T_start = TEvaIn_nominal, X_start = Medium.X_default = 0.01
u0 = [sys.sinSpeDX.vol.T_mix => 292.55, sys.sinSpeDX.vol.w_mix => 0.01 / 0.99]
prob = ODEProblem(sys, u0, (0.0, 86400.0); guesses = [sys.sinSpeDX.vol.p_mix => 101325.0])

tgrid, ref = read_mbl_reference(joinpath(DIR,
    "Buildings_Fluid_DXSystems_Cooling_AirSource_Validation_SingleSpeed.txt"))
tstops = sort(unique(vcat(kn1, kn2, kn3, kn4, kn5)))
sol = solve(prob, Rodas5P(); saveat = tgrid, tstops = tstops, reltol = 1e-6, abstol = 1e-8)
println("retcode = $(sol.retcode), points = $(length(sol.t))")

sim = (
    T = sol[sys.sinSpeDX.vol.T_mix] .- 273.15,
    X = sol[sys.sinSpeDX.vol.w_mix] ./ (1 .+ sol[sys.sinSpeDX.vol.w_mix]),
    m = sol[sys.sinSpeDX.preDro.mflow] .* (1 .+ sol[sys.sinSpeDX.preDro.w_in]),
    Q = sol[sys.sinSpeDX.Qflow.u],
    P = sol[sys.sinSpeDX.P.u],
)
mbl = (
    T = ref["sinSpeDX.vol.T"] .- 273.15,
    X = ref["sinSpeDX.vol.Xi[1]"],
    m = ref["sinSpeDX.m_flow"],
    Q = ref["QCoo_flow.y"],
    P = ref["sinSpeDX.P"],
)

rmse(a, b) = sqrt(mean((a .- b) .^ 2))
mbe(a, b) = mean(a .- b)
println("\nThis work vs MBL, whole day, $(length(tgrid)) points")
for (k, lbl) in ((:T, "T [degC]"), (:X, "X [kg/kg]"), (:m, "m [kg/s]"), (:Q, "Q [W]"), (:P, "P [W]"))
    @printf("  %-10s RMSE = %.3g   MBE = %.3g\n", lbl, rmse(sim[k], mbl[k]), mbe(sim[k], mbl[k]))
end

th = tgrid ./ 3600
function panel(k, ylab, unit, legpos)
    lbl = @sprintf("This work (RMSE: %.2g %s, MBE: %.2g %s)",
        rmse(sim[k], mbl[k]), unit, mbe(sim[k], mbl[k]), unit)
    p = plot(th, sim[k]; label = lbl, color = :green, linewidth = 3, ylabel = ylab,
        xlims = (0, 24), xticks = 0:4:24, legend = legpos,
        background_color_legend = Plots.RGBA(1, 1, 1, 0.6))
    plot!(p, th, mbl[k]; label = "MBL", color = :blue, linewidth = 1.5, linestyle = :dash)
    return p
end
fig = plot(
    panel(:T, "Outlet T [°C]", "°C", :bottomright),
    panel(:X, "Outlet X [kg/kg]", "kg/kg", :topright),
    panel(:Q, "Cooling [W]", "W", :topright),
    panel(:P, "Power [W]", "W", :bottomright),
    layout = (4, 1), size = (700, 900), xlabel = "Time [h]", left_margin = 5Plots.mm,
)
display(fig)
savefig(fig, joinpath(DIR, "result_images", "DXCoolingSingleSpeed.svg"))
