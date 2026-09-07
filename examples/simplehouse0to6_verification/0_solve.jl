using Plots
using Plots: mm, RGBA
default(legendfontsize = 10)
using Statistics
using Printf
using DifferentialEquations
using DataFrames
include(joinpath(@__DIR__, "helpers.jl"))
include(joinpath(@__DIR__, "0_model.jl"))

# MBL `SimpleHouse0` uses only `TDryBul`. `HGloHor` is checked here as well.

@mtkcompile sys = SimpleHouse0(df_weather = df_weather)

csv_path = joinpath(@__DIR__, "MBLresult_simplehouse0.csv")
isfile(csv_path) || error("Required CSV not found: $(csv_path)")

df_modelica = load_modelica_csv(csv_path)
modelica_time_s = df_modelica.time_s

tspan = (modelica_time_s[1], modelica_time_s[end])
prob = ODEProblem(sys, Pair[], tspan)
sol = solve(prob, Tsit5(); saveat = modelica_time_s, reltol = 1e-6)

TOut_sim_C = unit_K2C.(sol[sys.TOut.T.u])
HGlo_sim   = sol[sys.gaiHGlo.output.u]


# ── Load Modelica columns ──────────────────────────────────────────────────────
TOut_mbl_C  = df_modelica[!, "weaBus.TDryBul|degC"]
HGlo_mbl    = df_modelica[!, "weaBus.HGloHor|W/m2"]

mask       = modelica_time_s .< modelica_time_s[end]
TOut_sim_C = TOut_sim_C[mask]
TOut_mbl_C = TOut_mbl_C[mask]
HGlo_sim   = HGlo_sim[mask]
HGlo_mbl   = HGlo_mbl[mask]
time_hr    = modelica_time_s[mask] ./ 3600.0

# ── Metrics ───────────────────────────────────────────────────────────────────
rmse_TOut  = sqrt(mean((TOut_sim_C .- TOut_mbl_C) .^ 2))
mbe_TOut   = mean(TOut_sim_C .- TOut_mbl_C)
rmse_HGlo  = sqrt(mean((HGlo_sim  .- HGlo_mbl)   .^ 2))
mbe_HGlo   = mean(HGlo_sim  .- HGlo_mbl)

xtk = collect(range(0, 8760, length = 7))

# ── Plot: dry-bulb temperature ─────────────────────────────────────────────────
p_T_top = plot(
    time_hr, TOut_sim_C;
    label     = @sprintf("This work (RMSE: %.1e °C, MBE: %.1e °C)", rmse_TOut, mbe_TOut),
    title     = "Outdoor dry-bulb temperature (SimpleHouse0)",
    ylabel    = "Temperature [°C]",
    xlims     = (0, 8760),
    xticks    = (xtk, fill("", length(xtk))),
    color                    = :green,
    linewidth                = 3,
    legend                   = :bottom,
    background_color_legend  = RGBA(1, 1, 1, 0.6),
)
plot!(p_T_top, time_hr, TOut_mbl_C;
    label     = "MBL",
    color     = :blue,
    linewidth = 1,
    linestyle = :dash,
)

abserr_T = (TOut_sim_C .- TOut_mbl_C) .* 1e5   # shown in 1e-5 °C

p_T_bot = plot(
    time_hr, abserr_T;
    label   = "Error",
    ylabel  = "Error [10⁻⁵ °C]",
    xlabel  = "Time [hr]",
    xlims   = (0, 8760),
    xticks  = xtk,
    color                   = :gray,
    legend                  = :topleft,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
# hline!(p_T_bot, [0.0]; color = :black, linestyle = :dot, label = "")

fig_T = plot(p_T_top, p_T_bot, layout = grid(2, 1, heights = [0.67, 0.33]), size = (600, 550),
    left_margin = 0mm, right_margin = 3mm)
display(fig_T)

# ── Plot: global horizontal irradiance ────────────────────────────────────────
p_H_top = plot(
    time_hr, HGlo_sim;
    label     = @sprintf("This work (RMSE: %.1e W/m², MBE: %.1e W/m²)", rmse_HGlo, mbe_HGlo),
    title     = "Global horizontal irradiance (SimpleHouse0)",
    ylabel    = "Irradiance [W/m²]",
    xlims     = (0, 8760),
    xticks    = (xtk, fill("", length(xtk))),
    color                    = :green,
    linewidth                = 3,
    legend                   = :bottom,
    background_color_legend  = RGBA(1, 1, 1, 0.6),
)
plot!(p_H_top, time_hr, HGlo_mbl;
    label     = "MBL",
    color     = :blue,
    linewidth = 1,
    linestyle = :dash,
)

abserr_H = (HGlo_sim .- HGlo_mbl) .* 1e5   # shown in 1e-5 W/m²

p_H_bot = plot(
    time_hr, abserr_H;
    label   = "Error",
    ylabel  = "Error [10⁻⁵ W/m²]",
    xlabel  = "Time [hr]",
    xlims   = (0, 8760),
    xticks  = xtk,
    color                   = :gray,
    legend                  = :bottom,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
# hline!(p_H_bot, [0.0]; color = :black, linestyle = :dot, label = "")

fig_H = plot(p_H_top, p_H_bot, layout = grid(2, 1, heights = [0.67, 0.33]), size = (600, 550),
    left_margin = 0mm, right_margin = 3mm)
display(fig_H)
