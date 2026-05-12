using BuildEML
using DifferentialEquations
using Plots
using Plots: mm, RGBA
default(legendfontsize = 10)
using Statistics
using Printf
using DataFrames
include(joinpath(@__DIR__, "helpers.jl"))
include(joinpath(@__DIR__, "2_heatingonly_model.jl"))

epw_path = "src/Disturbances/Weather/weatherfile/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.epw"
df_weather = ReadEPW(epw_path)

@mtkbuild sys = SimpleHouseHeatingOnly(df_weather = df_weather)

@show unknowns(sys);

u0 = [
    sys.house.zon.T_mix => 293.15,
    sys.house.walCap.T => 293.15,
    sys.house.zon.w_mix => 0.008,
    sys.loop.rad.elems_1.T_mix => 293.15,
    sys.loop.rad.elems_2.T_mix => 293.15,
    sys.loop.rad.elems_3.T_mix => 293.15,
    sys.loop.rad.elems_4.T_mix => 293.15,
    sys.loop.rad.elems_5.T_mix => 293.15,
    sys.ctrller_hea.y => 1.0,
]

tspan = build_tspan(1, 1, 12, 31)
csv_path = joinpath(@__DIR__, "MBLresult_heatingonly.csv")
isfile(csv_path) || error("Required CSV not found: $(csv_path)")

df_modelica = load_modelica_csv(csv_path)
modelica_time_s = df_modelica.time_s
tzone_modelica_C = df_modelica[!, "zon.T|degC"]

prob = ODEProblem(sys, u0, tspan)
sol = solve(prob, Rodas5P(); saveat = modelica_time_s, abstol = 1e-6, reltol = 1e-4)

tzone_sim_C = unit_K2C.(sol[sys.house.RA_tap.T.u])
hea_signal  = sol[sys.ctrller_hea.y]
df_sim = DataFrame(
    time_s      = round.(Int, sol.t),
    Tzone_sim_C = tzone_sim_C,
    hea_on      = hea_signal .> 0.5,
)
df_mbl = DataFrame(
    time_s      = round.(Int, modelica_time_s),
    Tzone_mbl_C = tzone_modelica_C,
)
df_cmp = innerjoin(df_sim, df_mbl, on = :time_s)
nrow(df_cmp) == 0 && error("No matched timestamps between `sol.t` and Modelica CSV.")

rmse_tzone = sqrt(mean((df_cmp.Tzone_sim_C .- df_cmp.Tzone_mbl_C) .^ 2))
mbe_tzone  = mean(df_cmp.Tzone_sim_C .- df_cmp.Tzone_mbl_C)

t_plot    = df_cmp.time_s ./ 3600.0
t_end_plt = 8760.0
xtk       = collect(range(0, t_end_plt, length = 7))
xtk_labels = string.(round.(Int, xtk))

p_top = plot(
    t_plot, df_cmp.Tzone_sim_C;
    label  = @sprintf("This work (RMSE: %.2f °C, MBE: %.2f °C)", rmse_tzone, mbe_tzone),
    title  = "Zone air temperature (SimpleHouse, heating only)",
    ylabel = "Temperature [°C]",
    xlims  = (0, t_end_plt),
    xticks = (xtk, fill("", length(xtk))),
    color  = :green,
    linewidth = 3,
    legend = (0.4, 0.15),
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
plot!(p_top, t_plot, df_cmp.Tzone_mbl_C;
    label     = "MBL",
    color     = :blue,
    linewidth = 1,
    linestyle = :dash,
)

abserr = df_cmp.Tzone_sim_C .- df_cmp.Tzone_mbl_C

win_ma = 50
ma = [mean(abserr[max(1, i - win_ma ÷ 2):min(end, i + win_ma ÷ 2)]) for i in eachindex(abserr)]

p_bot = plot(
    t_plot, abserr;
    label  = "Pointwise error",
    ylabel = "Error [°C]",
    xlabel = "Time [hr]",
    xlims  = (0, t_end_plt),
    xticks = (xtk, xtk_labels),
    color  = :lightgray,
    linewidth = 1,
    legend = :topleft,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
plot!(p_bot, t_plot, ma;
    label     = "Moving average (+/-$(win_ma ÷ 2) pts)",
    color     = :red,
    linewidth = 2,
)
hline!(p_bot, [0.0]; color = :black, linestyle = :dot, linewidth = 1, label = "")

plot(p_top, p_bot, layout = grid(2, 1, heights = [0.67, 0.33]), size = (600, 550),
    left_margin = 0mm, right_margin = 3mm)
