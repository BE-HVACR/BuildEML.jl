using CSV, DataFrames

"""
    load_modelica_csv(csv_path) -> DataFrame

Load a Modelica result CSV, convert the time column to seconds if needed,
sort by time, and deduplicate by keeping the last entry per timestamp.
The first column is renamed to `time_s`; all other columns are preserved as-is.
"""
function load_modelica_csv(csv_path::String)
    df_raw = DataFrame(CSV.File(csv_path))
    time_col_name = lowercase(String(names(df_raw)[1]))
    time_raw = Float64.(df_raw[!, 1])
    time_s = occursin("|h", time_col_name) ? (time_raw .* 3600.0) : time_raw

    df = DataFrame(time_s = time_s)
    for col in names(df_raw)[2:end]
        df[!, col] = Float64.(df_raw[!, col])
    end

    sort!(df, :time_s)
    df = combine(groupby(df, :time_s),
        [col => last => col for col in names(df)[2:end]]...)
    return df
end
