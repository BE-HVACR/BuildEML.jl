using DifferentialEquations
using Plots
using Plots: mm, RGBA
default(legendfontsize = 10)
using Statistics
using Printf
using DataFrames
include(joinpath(@__DIR__, "helpers.jl"))
include(joinpath(@__DIR__, "3_model.jl"))

@mtkbuild sys = SimpleHouse3(df_weather = df_weather)
@show unknowns(sys)

u0 = [
    sys.walCap.T  => 293.15,
    sys.zon.T_mix => 293.15,
    sys.zon.w_mix => 0.008,
]

tspan = build_tspan(1, 1, 12, 31)
csv_path = joinpath(@__DIR__, "MBLresult_simplehouse3.csv")
isfile(csv_path) || error("Required CSV not found: $(csv_path)")

df_modelica = load_modelica_csv(csv_path)
modelica_time_s = df_modelica.time_s

prob = ODEProblem(sys, u0, tspan)
sol  = solve(prob, Rodas5P(); saveat = modelica_time_s, reltol = 1e-6)

walT_sim_C  = unit_K2C.(sol[sys.walCap.T])
zonT_sim_C  = unit_K2C.(sol[sys.RA_tap.T.u])
walT_mbl_C  = df_modelica[!, "walCap.T|degC"]
zonT_mbl_C  = df_modelica[!, "zon.T|degC"]

df_sim = DataFrame(
    time_s     = round.(Int, sol.t),
    walT_sim_C = walT_sim_C,
    zonT_sim_C = zonT_sim_C,
)
df_mbl = DataFrame(
    time_s     = round.(Int, modelica_time_s),
    walT_mbl_C = walT_mbl_C,
    zonT_mbl_C = zonT_mbl_C,
)
df_cmp = innerjoin(df_sim, df_mbl, on = :time_s)
nrow(df_cmp) == 0 && error("No matched timestamps between sol.t and Modelica CSV.")

mask   = df_cmp.time_s .< df_cmp.time_s[end]
df_cmp  = df_cmp[mask, :]
time_hr = df_cmp.time_s ./ 3600.0

# ── Metrics ───────────────────────────────────────────────────────────────────
rmse_walT = sqrt(mean((df_cmp.walT_sim_C .- df_cmp.walT_mbl_C) .^ 2))
mbe_walT  = mean(df_cmp.walT_sim_C .- df_cmp.walT_mbl_C)
rmse_zonT = sqrt(mean((df_cmp.zonT_sim_C .- df_cmp.zonT_mbl_C) .^ 2))
mbe_zonT  = mean(df_cmp.zonT_sim_C .- df_cmp.zonT_mbl_C)

xtk = collect(range(0, 8760, length = 7))

# ── Plot: wall temperature ─────────────────────────────────────────────────────
p_wal_top = plot(
    time_hr, df_cmp.walT_sim_C;
    label                    = @sprintf("This work (RMSE: %.2f°C, MBE: %.2f°C)", rmse_walT, mbe_walT),
    title                    = "Wall capacitor temperature (SimpleHouse3)",
    ylabel                   = "Temperature [°C]",
    xlims                    = (0, 8760),
    xticks                   = (xtk, fill("", length(xtk))),
    color                    = :green,
    linewidth                = 3,
    legend                   = :bottom,
    background_color_legend  = RGBA(1, 1, 1, 0.6),
)
plot!(p_wal_top, time_hr, df_cmp.walT_mbl_C;
    label     = "MBL",
    color     = :blue,
    linewidth = 1,
    linestyle = :dash,
)

relerr_wal = (df_cmp.walT_sim_C .- df_cmp.walT_mbl_C) ./ abs.(df_cmp.walT_mbl_C) .* 100

p_wal_bot = plot(
    time_hr, relerr_wal;
    label                   = "Relative error",
    ylabel                  = "Rel. Error [%]",
    xlabel                  = "Time [hr]",
    xlims                   = (0, 8760),
    ylims                   = (-1, 1),
    xticks                  = xtk,
    color                   = :gray,
    legend                  = :topleft,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)

plot(p_wal_top, p_wal_bot, layout = grid(2, 1, heights = [0.67, 0.33]), size = (600, 550),
    left_margin = 0mm, right_margin = 3mm)

# ── Plot: zone temperature ─────────────────────────────────────────────────────
p_zon_top = plot(
    time_hr, df_cmp.zonT_sim_C;
    label                    = @sprintf("This work (RMSE: %.2f°C, MBE: %.2f°C)", rmse_zonT, mbe_zonT),
    title                    = "Zone air temperature (SimpleHouse3)",
    ylabel                   = "Temperature [°C]",
    xlims                    = (0, 8760),
    xticks                   = (xtk, fill("", length(xtk))),
    color                    = :green,
    linewidth                = 3,
    legend                   = (0.4, 0.15),
    background_color_legend  = RGBA(1, 1, 1, 0.6),
)
plot!(p_zon_top, time_hr, df_cmp.zonT_mbl_C;
    label     = "MBL",
    color     = :blue,
    linewidth = 1,
    linestyle = :dash,
)

relerr_zon = (df_cmp.zonT_sim_C .- df_cmp.zonT_mbl_C) ./ abs.(df_cmp.zonT_mbl_C) .* 100

p_zon_bot = plot(
    time_hr, relerr_zon;
    label                   = "Relative error",
    ylabel                  = "Rel. Error [%]",
    xlabel                  = "Time [hr]",
    xlims                   = (0, 8760),
    ylims                   = (-1, 1),
    xticks                  = xtk,
    color                   = :gray,
    legend                  = :topleft,
    background_color_legend = RGBA(1, 1, 1, 0.6),
)

plot(p_zon_top, p_zon_bot, layout = grid(2, 1, heights = [0.67, 0.33]), size = (600, 550),
    left_margin = 0mm, right_margin = 3mm)
