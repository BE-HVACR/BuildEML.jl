using DifferentialEquations
using ModelingToolkit
using Plots
using Plots: mm, RGBA
default(legendfontsize = 10)
using Statistics
using Printf
using DataFrames
include(joinpath(@__DIR__, "helpers.jl"))
include(joinpath(@__DIR__, "6_model.jl"))

function sample_scalar(sol, sym, ts)
    [sol(t; idxs = sym) for t in ts]
end

function run_segmented_hysteresis(prob, sys, modelica_time_s; scan_dt = 60.0, y_off = 1e-4, y_on = 1.0)
    idx_yVent = findfirst(s -> isequal(s, sys.yVent), unknowns(sys))
    idx_zonT = findfirst(s -> isequal(s, sys.zon.T_mix), unknowns(sys))
    idx_yVent === nothing && error("Could not find sys.yVent in unknowns(sys).")
    idx_zonT === nothing && error("Could not find sys.zon.T_mix in unknowns(sys).")

    T_low_K = prob.ps[sys.T_low_vent]
    T_high_K = prob.ps[sys.T_high_vent]
    tf = prob.tspan[2]

    t_curr = prob.tspan[1]
    u_curr = copy(prob.u0)
    y_curr = y_off
    u_curr[idx_yVent] = y_curr

    seg_event_t = Float64[]
    seg_event_y = Float64[]
    out_t = Float64[]
    out_T = Float64[]
    out_y = Float64[]

    for _ in 1:100
        seg_out_times = modelica_time_s[(modelica_time_s .>= t_curr) .& (modelica_time_s .<= tf)]
        scan_times = collect(t_curr:scan_dt:tf)
        isempty(scan_times) && (scan_times = [t_curr])
        scan_times[end] == tf || push!(scan_times, tf)
        save_times = sort(unique(vcat(seg_out_times, scan_times)))

        seg_prob = remake(prob; u0 = copy(u_curr), tspan = (t_curr, tf))
        seg_sol = solve(seg_prob, Rodas5P(); saveat = save_times, tstops = save_times, reltol = 1e-6)

        T_scan_C = unit_K2C.(sample_scalar(seg_sol, sys.zon.T_mix, scan_times))
        hit_idx = if y_curr < 0.5
            findfirst(>=(25.0), T_scan_C)
        else
            findfirst(<=(23.0), T_scan_C)
        end
        t_event = isnothing(hit_idx) ? tf : scan_times[hit_idx]

        seg_keep_times = modelica_time_s[(modelica_time_s .>= t_curr) .& (modelica_time_s .< t_event)]
        append!(out_t, seg_keep_times)
        append!(out_T, unit_K2C.(sample_scalar(seg_sol, sys.RA_tap.T.u, seg_keep_times)))
        append!(out_y, fill(y_curr, length(seg_keep_times)))

        isnothing(hit_idx) && break

        push!(seg_event_t, t_event)
        y_curr = y_curr < 0.5 ? y_on : y_off
        push!(seg_event_y, y_curr)

        u_event = copy(seg_sol(t_event))
        u_event[idx_yVent] = y_curr
        u_curr = u_event
        t_curr = t_event
        t_curr >= tf && break
    end

    remaining = setdiff(modelica_time_s, out_t)
    if !isempty(remaining)
        final_prob = remake(prob; u0 = copy(u_curr), tspan = (t_curr, tf))
        final_sol = solve(final_prob, Rodas5P(); saveat = remaining, tstops = remaining, reltol = 1e-6)
        append!(out_t, remaining)
        append!(out_T, unit_K2C.(final_sol[sys.RA_tap.T.u]))
        append!(out_y, fill(y_curr, length(remaining)))
    end

    ord = sortperm(out_t)
    return (time_s = out_t[ord], zonT_C = out_T[ord], dam = out_y[ord], event_t = seg_event_t, event_y = seg_event_y)
end

@mtkbuild sys = SimpleHouse6(df_weather = df_weather)

u0 = [
    sys.walCap.T           => 293.15,
    sys.zon.T_mix          => 293.15,
    sys.zon.w_mix          => 0.008,
    sys.rad.elems_1.T_mix  => 333.15,
    sys.rad.elems_2.T_mix  => 328.15,
    sys.rad.elems_3.T_mix  => 323.15,
    sys.rad.elems_4.T_mix  => 318.15,
    sys.rad.elems_5.T_mix  => 313.15,
    sys.ctrllerHea.y       => 1.0,
    sys.yVent              => 1e-4,
]

tspan = (0.0, 1e6)
csv_path = joinpath(@__DIR__, "MBLresult_simplehouse6.csv")
isfile(csv_path) || error("Required CSV not found: $(csv_path)")

df_modelica = load_modelica_csv(csv_path)
modelica_time_s = df_modelica.time_s
guesses = [
    sys.hexRec.m1_flow => 0.0,
    sys.hexRec.m2_flow => 0.0,
    sys.hexRec.T1_out  => 293.15,
    sys.hexRec.T2_out  => 293.15,
]
prob = ODEProblem(sys, u0, tspan; guesses = guesses)

sim = run_segmented_hysteresis(prob, sys, modelica_time_s; scan_dt = 60.0)

zonT_sim_C = sim.zonT_C
dam_sim = sim.dam
zonT_mbl_C = df_modelica[!, "zon.T|degC"]
dam_mbl = df_modelica[!, "vavDam.y"]

df_sim = DataFrame(
    time_s     = round.(Int, sim.time_s),
    zonT_sim_C = zonT_sim_C,
    dam_sim    = dam_sim,
)
df_mbl = DataFrame(
    time_s     = round.(Int, modelica_time_s),
    zonT_mbl_C = zonT_mbl_C,
    dam_mbl    = dam_mbl,
)
df_cmp = innerjoin(df_sim, df_mbl, on = :time_s)
nrow(df_cmp) == 0 && error("No matched timestamps between simulation and Modelica CSV.")

time_s = Float64.(df_cmp.time_s)
rmse_zonT = sqrt(mean((df_cmp.zonT_sim_C .- df_cmp.zonT_mbl_C) .^ 2))
mbe_zonT  = mean(df_cmp.zonT_sim_C .- df_cmp.zonT_mbl_C)

t_max     = time_s[end]
use_hours = t_max > 1e6

if use_hours
    t_plot     = time_s ./ 3600.0
    t_end_plt  = t_max  / 3600.0
    xtk        = collect(range(0, t_end_plt, length = 7))
    xtk_labels = string.(round.(Int, xtk))
    xlabel_str = "Time [hr]"
else
    t_plot     = time_s
    t_end_plt  = t_max
    xtk        = collect(range(0, t_end_plt, length = 6))
    xtk_labels = ["0", "2×10⁵", "4×10⁵", "6×10⁵", "8×10⁵", "1×10⁶"]
    xlabel_str = "Time [s]"
end

p_top = plot(
    t_plot, df_cmp.zonT_sim_C;
    label                   = @sprintf("This work (RMSE: %.2f°C, MBE: %.2f°C)", rmse_zonT, mbe_zonT),
    title                   = "Zone air temperature (SimpleHouse6)",
    ylabel                  = "Temperature [°C]",
    xlims                   = (0, t_end_plt),
    xticks                  = (xtk, fill("", length(xtk))),
    color                   = :green,
    linewidth               = 3,
    legend                  = (0.4, 0.15),
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
plot!(p_top, t_plot, df_cmp.zonT_mbl_C;
    label     = "MBL",
    color     = :blue,
    linewidth = 1,
    linestyle = :dash,
)

abserr = df_cmp.zonT_sim_C .- df_cmp.zonT_mbl_C
win_ma = 50
ma = [mean(abserr[max(1, i - win_ma ÷ 2):min(end, i + win_ma ÷ 2)]) for i in eachindex(abserr)]

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
    label     = "Moving average (±$(win_ma ÷ 2) pts)",
    color     = :red,
    linewidth = 2,
)
hline!(p_bot, [0.0]; color = :black, linestyle = :dot, linewidth = 1, label = "")

display(plot(p_top, p_bot, layout = grid(2, 1, heights = [0.67, 0.33]), size = (600, 550),
    left_margin = 0mm, right_margin = 3mm))

rmse_dam = sqrt(mean((df_cmp.dam_sim .- df_cmp.dam_mbl) .^ 2))
mbe_dam  = mean(df_cmp.dam_sim .- df_cmp.dam_mbl)

p_dam_top = plot(
    t_plot, df_cmp.dam_sim;
    label                   = @sprintf("This work (RMSE: %.2f, MBE: %.2f)", rmse_dam, mbe_dam),
    title                   = "Damper opening signal (SimpleHouse6)",
    ylabel                  = "Opening [-]",
    xlims                   = (0, t_end_plt),
    xticks                  = (xtk, fill("", length(xtk))),
    color                   = :green,
    linewidth               = 3,
    legend                  = :topleft,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)
plot!(p_dam_top, t_plot, df_cmp.dam_mbl;
    label     = "Modelica",
    color     = :blue,
    linewidth = 1,
    linestyle = :dash,
)

abserr_dam = df_cmp.dam_sim .- df_cmp.dam_mbl
p_dam_bot = plot(
    t_plot, abserr_dam;
    label                   = "Absolute error",
    ylabel                  = "Error [-]",
    xlabel                  = xlabel_str,
    xlims                   = (0, t_end_plt),
    xticks                  = (xtk, xtk_labels),
    color                   = :gray,
    legend                  = :topleft,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)

plot(p_dam_top, p_dam_bot, layout = grid(2, 1, heights = [0.67, 0.33]), size = (600, 550),
    left_margin = 0mm, right_margin = 3mm)

