using DifferentialEquations
using Plots
using Plots: mm, RGBA
default(legendfontsize = 10)
using Statistics
using Printf
using DataFrames
include(joinpath(@__DIR__, "helpers.jl"))
include(joinpath(@__DIR__, "5_model.jl"))

@mtkbuild sys = SimpleHouse5(df_weather = df_weather)
# @show unknowns(sys)

u0 = [
    sys.walCap.T          => 293.15,
    sys.zon.T_mix         => 293.15,
    sys.zon.w_mix         => 0.008,
    sys.rad.elems_1.T_mix => 333.15,
    sys.rad.elems_2.T_mix => 328.15,
    sys.rad.elems_3.T_mix => 323.15,
    sys.rad.elems_4.T_mix => 318.15,
    sys.rad.elems_5.T_mix => 313.15,
    sys.ctrllerHea.y      => 1.0,
]

tspan = (0.0, 1e6) 
# tspan = build_tspan(1, 1, 12, 31)
csv_path = joinpath(@__DIR__, "MBLresult_simplehouse5.csv")
isfile(csv_path) || error("Required CSV not found: $(csv_path)")

df_modelica = load_modelica_csv(csv_path)
modelica_time_s = df_modelica.time_s

prob = ODEProblem(sys, u0, tspan)
sol  = solve(prob, Rodas5P(); saveat = modelica_time_s, reltol = 1e-6)

zonT_sim_C = unit_K2C.(sol[sys.RA_tap.T.u])
zonT_mbl_C = df_modelica[!, "zon.T|degC"]

df_sim = DataFrame(time_s = round.(Int, sol.t),           zonT_sim_C = zonT_sim_C)
df_mbl = DataFrame(time_s = round.(Int, modelica_time_s), zonT_mbl_C = zonT_mbl_C)
df_cmp = innerjoin(df_sim, df_mbl, on = :time_s)
nrow(df_cmp) == 0 && error("No matched timestamps between sol.t and Modelica CSV.")

time_s = Float64.(df_cmp.time_s)

# ── Metrics ───────────────────────────────────────────────────────────────────
rmse_zonT = sqrt(mean((df_cmp.zonT_sim_C .- df_cmp.zonT_mbl_C) .^ 2))
mbe_zonT  = mean(df_cmp.zonT_sim_C .- df_cmp.zonT_mbl_C)

t_max     = time_s[end]
use_hours = t_max > 1e6   # full-year run → hours axis

if use_hours
    t_plot    = time_s ./ 3600.0
    t_end_plt = t_max  / 3600.0
    xtk       = collect(range(0, t_end_plt, length = 7))
    xtk_labels = string.(round.(Int, xtk))
    xlabel_str = "Time [hr]"
else
    t_plot    = time_s
    t_end_plt = t_max
    xtk       = collect(range(0, t_end_plt, length = 6))
    xtk_labels = ["0", "2×10⁵", "4×10⁵", "6×10⁵", "8×10⁵", "1×10⁶"]
    xlabel_str = "Time [s]"
end

# ── Plot ──────────────────────────────────────────────────────────────────────
p_top = plot(
    t_plot, df_cmp.zonT_sim_C;
    label                    = @sprintf("This work (RMSE: %.2f°C, MBE: %.2f°C)", rmse_zonT, mbe_zonT),
    title                    = "Zone air temperature (SimpleHouse5)",
    ylabel                   = "Temperature [°C]",
    xlims                    = (0, t_end_plt),
    xticks                   = (xtk, fill("", length(xtk))),
    color                    = :green,
    linewidth                = 3,
    legend                   = (0.4, 0.15),
    background_color_legend  = RGBA(1, 1, 1, 0.6),
)
plot!(p_top, t_plot, df_cmp.zonT_mbl_C;
    label     = "MBL",
    color     = :blue,
    linewidth = 1,
    linestyle = :dash,
)

abserr = df_cmp.zonT_sim_C .- df_cmp.zonT_mbl_C

# Moving average over a window of ~50 points
win_ma = 50
ma = [mean(abserr[max(1,i-win_ma÷2):min(end,i+win_ma÷2)]) for i in eachindex(abserr)]

p_bot = plot(
    t_plot, abserr;
    label                   = "Pointwise error",
    ylabel                  = "Error [°C]",
    xlabel                  = xlabel_str,
    xlims                   = (0, t_end_plt),
    xticks                  = (xtk, xtk_labels),
    color                   = :lightgray,
    linewidth               = 1,
    legend                  = :topleft,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
plot!(p_bot, t_plot, ma;
    label     = "Moving average (±$(win_ma÷2) pts)",
    color     = :red,
    linewidth = 2,
)
hline!(p_bot, [0.0]; color = :black, linestyle = :dot, linewidth = 1, label = "")

plot(p_top, p_bot, layout = grid(2, 1, heights = [0.67, 0.33]), size = (600, 550),
    left_margin = 0mm, right_margin = 3mm)
