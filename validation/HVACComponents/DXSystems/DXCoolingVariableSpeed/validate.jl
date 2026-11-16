# Validation of DXCoolingVariableSpeed against MBL v11
# Buildings.Fluid.DXSystems.Cooling.AirSource.Examples.VariableSpeed.
#
# VariableSpeed.mo is the MBL example; its coil data and inputs are reproduced
# below and the speed-ratio table is parsed from it.
# Buildings_Fluid_DXSystems_Cooling_AirSource_Examples_VariableSpeed.txt is the MBL
# Dymola reference result file (Resources/ReferenceResults/Dymola), unchanged.
#
# Differences from the MBL model:
# - accumulated-water reevaporation is not modeled (MBL computes it, but the coil
#   does not switch off once it is on in this case)
# - the minimum-speed hysteresis relaxes to its target with a 10 ms time constant

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

src = read(joinpath(DIR, "VariableSpeed.mo"), String)
tSR, vSR = read_timetable(src, "speRat")
const ITP_SR = LinearInterpolation(vSR, tSR; extrapolation = ExtrapolationType.Constant)

in_SR(x) = ITP_SR(x)
# Ramp p: offset 101325, height dp_nominal = 1000, startTime 100, duration 600
in_P(x) = 101325.0 + 1000.0 * clamp((x - 100.0) / 600.0, 0.0, 1.0)
# Ramp TEvaIn: offset 293.15, height 5, startTime 900, duration 600
in_TE(x) = 293.15 + 5.0 * clamp((x - 900.0) / 600.0, 0.0, 1.0)
@register_symbolic in_SR(x)
@register_symbolic in_P(x)
@register_symbolic in_TE(x)

# datCoi of the MBL example; X of the source is Medium.X_default = 0.01
const STAGES = [
    dx_stage(spe = 900 / 60, Qflow_nominal = -12000.0, COP_nominal = 3.0,
        SHR_nominal = 0.8, mflow_nominal = 0.9, perf_curves = DX_CURVE_I),
    dx_stage(spe = 1200 / 60, Qflow_nominal = -18000.0, COP_nominal = 3.0,
        SHR_nominal = 0.8, mflow_nominal = 1.2, perf_curves = DX_CURVE_I),
    dx_stage(spe = 1800 / 60, Qflow_nominal = -21000.0, COP_nominal = 3.0,
        SHR_nominal = 0.8, mflow_nominal = 1.5, perf_curves = DX_CURVE_II),
    dx_stage(spe = 2400 / 60, Qflow_nominal = -30000.0, COP_nominal = 3.0,
        SHR_nominal = 0.8, mflow_nominal = 1.8, perf_curves = DX_CURVE_III),
]

@mtkmodel ValidationSource begin
    @components begin
        port = AirPort()
    end
    @equations begin
        port.p ~ in_P(t)
        port.T_ofo ~ in_TE(t)
        port.w_ofo ~ 0.01 / 0.99
    end
end

@mtkmodel VariableSpeedValidation begin
    @components begin
        sou = ValidationSource()
        varSpeDX = DXCoolingVariableSpeed(stages = STAGES, dp_nominal = 1000.0,
            if_hysteresis_shutoff = true)
        sin = AirBoundary(role = :sink, mode = :pressure)
    end
    @equations begin
        connect(sou.port, varSpeDX.port_a)
        connect(varSpeDX.port_b, sin.port)
        sin.p.u ~ 101325.0
        varSpeDX.speRat.u ~ in_SR(t)
        varSpeDX.TConIn.u ~ 298.15
    end
end

@mtkcompile sys = VariableSpeedValidation()

# MBL: T_start = TEvaIn_nominal of stage 1, X_start = Medium.X_default = 0.01, coil off
u0 = [sys.varSpeDX.vol.T_mix => 292.55, sys.varSpeDX.vol.w_mix => 0.01 / 0.99,
    sys.varSpeDX.on_sta => 0.0]
prob = ODEProblem(sys, u0, (0.0, 3600.0); guesses = [sys.varSpeDX.vol.p_mix => 101325.0])

tgrid, ref = read_mbl_reference(joinpath(DIR,
    "Buildings_Fluid_DXSystems_Cooling_AirSource_Examples_VariableSpeed.txt"))
tstops = sort(unique(vcat(tSR, [100.0, 700.0, 900.0, 1500.0])))
sol = solve(prob, Rodas5P(); saveat = tgrid, tstops = tstops, reltol = 1e-6, abstol = 1e-8)
println("retcode = $(sol.retcode), points = $(length(sol.t))")

sim = (
    T = sol[sys.varSpeDX.vol.T_mix] .- 273.15,
    X = sol[sys.varSpeDX.vol.w_mix] ./ (1 .+ sol[sys.varSpeDX.vol.w_mix]),
    Q = sol[sys.varSpeDX.core.Qflow.u],
    mWat = sol[sys.varSpeDX.core.mWat_flow.u],
    P = sol[sys.varSpeDX.P.u],
)
mbl = (
    T = ref["varSpeDX.vol.T"] .- 273.15,
    X = ref["varSpeDX.vol.X_w"],
    Q = ref["varSpeDX.dxCoi.Q_flow"],
    mWat = ref["varSpeDX.dxCoi.mWat_flow"],
    P = ref["varSpeDX.P"],
)

rmse(a, b) = sqrt(mean((a .- b) .^ 2))
mbe(a, b) = mean(a .- b)
println("\nThis work vs MBL, $(length(tgrid)) points")
for (k, lbl) in ((:T, "T [degC]"), (:X, "X [kg/kg]"), (:Q, "Q [W]"), (:mWat, "mWat [kg/s]"), (:P, "P [W]"))
    @printf("  %-12s RMSE = %.3g   MBE = %.3g   max|d| = %.3g\n", lbl,
        rmse(sim[k], mbl[k]), mbe(sim[k], mbl[k]), maximum(abs.(sim[k] .- mbl[k])))
end

tm = tgrid ./ 60
function panel(k, ylab, unit, legpos)
    lbl = @sprintf("This work (RMSE: %.2g %s, MBE: %.2g %s)",
        rmse(sim[k], mbl[k]), unit, mbe(sim[k], mbl[k]), unit)
    p = plot(tm, sim[k]; label = lbl, color = :green, linewidth = 3, ylabel = ylab,
        xlims = (0, 60), xticks = 0:10:60, legend = legpos,
        background_color_legend = Plots.RGBA(1, 1, 1, 0.6))
    plot!(p, tm, mbl[k]; label = "MBL", color = :blue, linewidth = 1.5, linestyle = :dash)
    return p
end
fig = plot(
    panel(:T, "Outlet T [°C]", "°C", :topright),
    panel(:X, "Outlet X [kg/kg]", "kg/kg", :topright),
    panel(:Q, "Cooling [W]", "W", :topright),
    panel(:P, "Power [W]", "W", :bottomright),
    layout = (4, 1), size = (700, 900), xlabel = "Time [min]", left_margin = 5Plots.mm,
)
display(fig)
savefig(fig, joinpath(DIR, "result_images", "DXCoolingVariableSpeed.svg"))
