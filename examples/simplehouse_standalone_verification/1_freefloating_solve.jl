using BuildEML
using DifferentialEquations
using Plots
using Plots: mm, RGBA
default(legendfontsize = 10)
using Statistics
using Printf
using DataFrames
include(joinpath(@__DIR__, "helpers.jl"))

epw_path = "src/Disturbances/Weather/weatherfile/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.epw"
df_weather = ReadEPW(epw_path)

@mtkbuild sys = SimpleHouseCore(
    df_weather = df_weather,
    n_airports = 2,
    port_role = [:inlet, :outlet],
    mAir_flow_infil = 1e-6,
    use_w_in = false
)

@show unknowns(sys);

u0 = [
    sys.zon.T_mix => 293.15,
    sys.walCap.T => 293.15,
    sys.zon.w_mix => 0.008,
]

tspan = build_tspan(1, 1, 12, 31)
csv_path = joinpath(@__DIR__, "MBLresult_freefloating.csv")
isfile(csv_path) || error("Required CSV not found: $(csv_path)")

df_modelica = load_modelica_csv(csv_path)
modelica_time_s  = df_modelica.time_s
Tzone_modelica_C = df_modelica[!, "zon.T|degC"]
Tout_modelica_C  = df_modelica[!, "TOut.port.T|degC"]

prob = ODEProblem(sys, u0, tspan)
sol = solve(prob, Rodas5P(); saveat = modelica_time_s, abstol = 1e-6, reltol = 1e-6)

Tzone_sim_C = unit_K2C.(sol[sys.RA_tap.T.u])
Tout_sim_C  = unit_K2C.(sol[sys.weaBus.TDryBul.u])

df_sim = DataFrame(
    time_s      = round.(Int, sol.t),
    Tzone_sim_C = Tzone_sim_C,
    Tout_sim_C  = Tout_sim_C,
)
df_mbl = DataFrame(
    time_s      = round.(Int, modelica_time_s),
    Tzone_mbl_C = Tzone_modelica_C,
    Tout_mbl_C  = Tout_modelica_C,
)
df_cmp = innerjoin(df_sim, df_mbl, on = :time_s)
nrow(df_cmp) == 0 && error("No matched timestamps between `sol.t` and Modelica CSV.")

mask   = df_cmp.time_s .< df_cmp.time_s[end]
df_cmp = df_cmp[mask, :]

rmse_Tzone = sqrt(mean((df_cmp.Tzone_sim_C .- df_cmp.Tzone_mbl_C) .^ 2))
mbe_Tzone  = mean(df_cmp.Tzone_sim_C .- df_cmp.Tzone_mbl_C)
rmse_Tout  = sqrt(mean((df_cmp.Tout_sim_C  .- df_cmp.Tout_mbl_C)  .^ 2))
mbe_Tout   = mean(df_cmp.Tout_sim_C  .- df_cmp.Tout_mbl_C)

t_plot     = df_cmp.time_s ./ 3600.0
t_end_plt  = 8760.0
xtk        = collect(range(0, t_end_plt, length = 7))
xtk_labels = string.(round.(Int, xtk))

p_top = plot(
    t_plot, df_cmp.Tzone_sim_C;
    label  = @sprintf("Tzone, this work (RMSE: %.2f °C, MBE: %.2f °C)", rmse_Tzone, mbe_Tzone),
    title  = "Zone air temperature (SimpleHouse, free-floating)",
    ylabel = "Temperature [°C]",
    xlims  = (0, t_end_plt),
    xticks = (xtk, fill("", length(xtk))),
    color  = :green,
    linewidth = 2,
    legend = (0.25, 0.28),
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
plot!(p_top, t_plot, df_cmp.Tzone_mbl_C;
    label     = "Tzone, MBL",
    color     = :blue,
    linewidth = 1,
    linestyle = :dash,
)
plot!(p_top, t_plot, df_cmp.Tout_sim_C;
    label     = @sprintf("Tout, this work (RMSE: %.2f °C, MBE: %.2f °C)", rmse_Tout, mbe_Tout),
    color     = :orange,
    linewidth = 2,
)
plot!(p_top, t_plot, df_cmp.Tout_mbl_C;
    label     = "Tout, MBL",
    color     = :purple,
    linewidth = 1,
    linestyle = :dash,
)

abserr_Tout  = df_cmp.Tout_sim_C  .- df_cmp.Tout_mbl_C
abserr_Tzone = df_cmp.Tzone_sim_C .- df_cmp.Tzone_mbl_C

p_mid = plot(
    t_plot, abserr_Tout;
    label  = @sprintf("Tout error (RMSE: %.2f °C, MBE: %.2f °C)", rmse_Tout, mbe_Tout),
    title  = "Tout error",
    ylabel = "Error [°C]",
    xlims  = (0, t_end_plt),
    xticks = (xtk, fill("", length(xtk))),
    color  = :orange,
    linewidth = 1,
    legend = :topleft,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)

p_bot = plot(
    t_plot, abserr_Tzone;
    label  = @sprintf("Tzone error (RMSE: %.2f °C, MBE: %.2f °C)", rmse_Tzone, mbe_Tzone),
    title  = "Tzone error",
    ylabel = "Error [°C]",
    xlabel = "Time [hr]",
    xlims  = (0, t_end_plt),
    xticks = (xtk, xtk_labels),
    ylims  = (-0.1, 0.1),
    color  = :green,
    linewidth = 1,
    legend = :topleft,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)

plot(p_top, p_mid, p_bot, layout = grid(3, 1, heights = [0.5, 0.25, 0.25]), size = (600, 750),
    left_margin = 5mm, right_margin = 3mm)
