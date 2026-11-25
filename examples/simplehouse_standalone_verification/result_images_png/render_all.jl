ENV["GKSwstype"] = "100"

project_dir = normpath(joinpath(@__DIR__, "..", "..", ".."))
script_path = @__FILE__

if isempty(ARGS)
    for case_name in ("1_freefloating", "2_heatingonly", "3_fulloperation")
        run(`$(Base.julia_cmd()) --project=$project_dir $script_path $case_name`)
    end
    exit()
end

case_name = only(ARGS)
solve_path = normpath(joinpath(@__DIR__, "..", "$(case_name)_solve.jl"))
include(solve_path)

svg_dir = mkpath(joinpath(@__DIR__, "..", "result_images_svg"))
function save_both(p, name)
    savefig(p, joinpath(svg_dir, "$(name).svg"))
    savefig(plot!(p; dpi = 300), joinpath(@__DIR__, "$(name).png"))
end

if case_name == "1_freefloating"
    save_both(fig, "1_standalone_freefloating")
elseif case_name == "2_heatingonly"
    save_both(fig, "2_standalone_heatingonly")
    save_both(p_zoom, "2_standalone_heatingonly_zoom12h")
elseif case_name == "3_fulloperation"
    save_both(p_cmp, "3_standalone_fulloperation")
else
    error("Unknown case: $case_name")
end
