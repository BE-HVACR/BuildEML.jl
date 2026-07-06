using BuildEML
using DifferentialEquations
using Plots
using Plots: mm, RGBA
default(legendfontsize = 10)
using Statistics
using Printf
using DataFrames
include(joinpath(@__DIR__, "helpers.jl"))
include(joinpath(@__DIR__, "3_fulloperation_model.jl"))

epw_path = "src/Disturbances/Weather/weatherfile/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.epw"
df_weather = ReadEPW(epw_path)

@mtkbuild sys = SimpleHouseTotal(df_weather = df_weather)

u0 = [
    sys.house.zon.T_mix => 293.15,
    sys.house.walCap.T => 293.15,
    sys.house.zon.w_mix => 0.008,
    sys.loop_hea.rad.elems_1.T_mix => 293.15,
    sys.loop_hea.rad.elems_2.T_mix => 293.15,
    sys.loop_hea.rad.elems_3.T_mix => 293.15,
    sys.loop_hea.rad.elems_4.T_mix => 293.15,
    sys.loop_hea.rad.elems_5.T_mix => 293.15,
    sys.ctrller_hea.y => 1.0,
]

guesses = [
    sys.loop_cool.dam.dp => 10.0,
    sys.loop_cool.dam.m_flow => 0.01,
    sys.loop_cool.dam.sgn_dp => 1.0,
    sys.loop_cool.dam.f_y => 0.5,
]

tspan = build_tspan(1, 1, 12, 31)
csv_path = joinpath(@__DIR__, "MBLresult_fulloperation.csv")
isfile(csv_path) || error("Required CSV not found: $(csv_path)")

df_modelica = load_modelica_csv(csv_path)
modelica_time_s = df_modelica.time_s
Tzone_modelica_C = df_modelica[!, "zon.T|degC"]

prob = ODEProblem(sys, u0, tspan; guesses = guesses)
sol = solve(prob, Rodas5P(); saveat = modelica_time_s, abstol = 1e-6, reltol = 1e-4)

Tzone_sim_C = unit_K2C.(sol[sys.house.RA_tap.T.u])

df_sim = DataFrame(
    time_s = round.(Int, sol.t),
    Tzone_sim_C = Tzone_sim_C,
)
df_mbl = DataFrame(
    time_s = round.(Int, modelica_time_s),
    Tzone_mbl_C = Tzone_modelica_C,
)
df_cmp = innerjoin(df_sim, df_mbl, on = :time_s)
nrow(df_cmp) == 0 && error("No matched timestamps between sol.t and Modelica CSV.")

time_s = Float64.(df_cmp.time_s)
rmse_Tzone = sqrt(mean((df_cmp.Tzone_sim_C .- df_cmp.Tzone_mbl_C) .^ 2))
mbe_Tzone  = mean(df_cmp.Tzone_sim_C .- df_cmp.Tzone_mbl_C)

t_max = time_s[end]
use_hours = t_max > 1e6

if use_hours
    t_plot = time_s ./ 3600.0
    t_end_plt = t_max / 3600.0
    xtk = collect(range(0, t_end_plt, length = 7))
    xtk_labels = string.(round.(Int, xtk))
    xlabel_str = "Time [hr]"
else
    t_plot = time_s
    t_end_plt = t_max
    xtk = collect(range(0, t_end_plt, length = 6))
    xtk_labels = ["0", "2e5", "4e5", "6e5", "8e5", "1e6"]
    xlabel_str = "Time [s]"
end

p_top = plot(
    t_plot, df_cmp.Tzone_sim_C;
    label = @sprintf("This work (RMSE: %.2f °C, MBE: %.2f °C)", rmse_Tzone, mbe_Tzone),
    title = "SimpleHouse, complete", # "Zone air temperature (SimpleHouse)",
    ylabel = "Temperature [°C]",
    xlims = (0, t_end_plt),
    xticks = (xtk, fill("", length(xtk))),
    color = :green,
    linewidth = 3,
    legend = (0.37, 0.15),
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
plot!(p_top, t_plot, df_cmp.Tzone_mbl_C;
    label = "MBL",
    color = :blue,
    linewidth = 1,
    linestyle = :dash,
)

abserr = df_cmp.Tzone_sim_C .- df_cmp.Tzone_mbl_C

win_ma = 50
ma = [mean(abserr[max(1, i - win_ma ÷ 2):min(end, i + win_ma ÷ 2)]) for i in eachindex(abserr)]

p_bot = plot(
    t_plot, abserr;
    label = "Pointwise error",
    ylabel = "Error [°C]",
    xlabel = xlabel_str,
    xlims = (0, t_end_plt),
    xticks = (xtk, xtk_labels),
    color = :lightgray,
    linewidth = 1,
    legend = :topleft,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
plot!(p_bot, t_plot, ma;
    label = "Moving average (+/-$(win_ma ÷ 2) pts)",
    color = :red,
    linewidth = 2,
)
hline!(p_bot, [0.0]; color = :black, linestyle = :dot, linewidth = 1, label = "")

#=
# The MBL SimpleHouse case does not expose heating/cooling source power for comparison.
# This diagnostic plot shows only H/C source load, power, PLR, and control, excluding fan and pump energy.
hvac_plot_mode = :all  # :all, :heating, or :cooling
@assert hvac_plot_mode in (:all, :heating, :cooling) "hvac_plot_mode must be :all, :heating, or :cooling."

p_cmp = plot(p_top, p_bot, layout = grid(2, 1, heights = [0.67, 0.33]), size = (600, 550),
    left_margin = 0mm, right_margin = 3mm)
display(p_cmp)

t_plot_sol = use_hours ? Float64.(sol.t) ./ 3600.0 : Float64.(sol.t)
heating_power_W = sol[sys.loop_hea.heatingsrc_power.power.u]
cooling_power_W = sol[sys.loop_cool.coolingsource.power.u]
heating_load_W = sol[sys.loop_hea.hea.Qflow.u]
cooling_load_W = sol[sys.loop_cool.coolingsource.Q_thermal.u]
heating_capacity_W = 700.0
heating_load_abs_W = max.(heating_load_W, 0.0)
cooling_load_abs_W = max.(-cooling_load_W, 0.0)
heating_plr = heating_load_abs_W ./ heating_capacity_W
cooling_plr = sol[sys.loop_cool.coolingsource.yPL.u]
heating_control = sol[sys.ctrller_hea.y]
cooling_control = sol[sys.loop_cool.conDam.ctr_output.u]
plot_heating = hvac_plot_mode in (:all, :heating)
plot_cooling = hvac_plot_mode in (:all, :cooling)

p_load = plot(;
    title = "H/C thermal load",
    ylabel = "Thermal load [W]",
    xlims = (0, t_end_plt),
    xticks = (xtk, fill("", length(xtk))),
    legend = :topright,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
if plot_heating
    plot!(p_load, t_plot_sol, heating_load_abs_W;
        label = "Heating",
        color = :red,
        linewidth = 2,
    )
end
if plot_cooling
    plot!(p_load, t_plot_sol, cooling_load_abs_W;
        label = "Cooling",
        color = :blue,
        linewidth = 2,
    )
end

p_power = plot(;
    title = "H/C electric power",
    ylabel = "Power [W]",
    xlims = (0, t_end_plt),
    xticks = (xtk, fill("", length(xtk))),
    legend = :topright,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
if plot_heating
    plot!(p_power, t_plot_sol, heating_power_W;
        label = "Heating",
        color = :red,
        linewidth = 2,
    )
end
if plot_cooling
    plot!(p_power, t_plot_sol, cooling_power_W;
        label = "Cooling",
        color = :blue,
        linewidth = 2,
    )
end

p_plr = plot(;
    title = "H/C part-load ratio",
    ylabel = "PLR [-]",
    xlabel = xlabel_str,
    xlims = (0, t_end_plt),
    xticks = (xtk, xtk_labels),
    legend = :topright,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
if plot_heating
    plot!(p_plr, t_plot_sol, heating_plr;
        label = "Heating",
        color = :red,
        linewidth = 2,
    )
end
if plot_cooling
    plot!(p_plr, t_plot_sol, cooling_plr;
        label = "Cooling",
        color = :blue,
        linewidth = 2,
    )
end
hline!(p_plr, [0.0]; color = :black, linestyle = :dot, linewidth = 1, label = "")

p_control = plot(;
    title = "H/C control",
    ylabel = "Control [-]",
    xlabel = xlabel_str,
    xlims = (0, t_end_plt),
    xticks = (xtk, xtk_labels),
    legend = :topright,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
if plot_heating
    plot!(p_control, t_plot_sol, heating_control;
        label = "Heating",
        color = :red,
        linewidth = 2,
    )
end
if plot_cooling
    plot!(p_control, t_plot_sol, cooling_control;
        label = "Cooling",
        color = :blue,
        linewidth = 2,
    )
end

p_hvac = plot(p_load, p_power, p_plr, p_control, layout = grid(2, 2), size = (1000, 650),
    left_margin = 3mm, right_margin = 3mm)
display(p_hvac)
=#