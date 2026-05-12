function _validate_schedule_dataframe(df::AbstractDataFrame; time_col::Symbol = :time)
    nameset = Set(Symbol.(names(df)))
    time_col in nameset || throw(ArgumentError("Schedule DataFrame is missing required time column $(time_col)."))

    missing_cols = Symbol[col for col in SCHEDULE_REQUIRED_COLUMNS if !(col in nameset)]
    isempty(missing_cols) || throw(ArgumentError("Schedule DataFrame is missing required columns: $(join(string.(missing_cols), ", "))."))

    nrow(df) > 0 || throw(ArgumentError("Schedule DataFrame must contain at least one row."))
    return df
end


"""
    ReadScheduleCSV(path; time_col=:time, default_step_seconds=3600.0)

Read a schedule CSV and validate that it contains the columns required by
`ScheduleBus`. If the `time` column is missing, it is generated assuming a
uniform step of `default_step_seconds`.
"""
function ReadScheduleCSV(csv_path::AbstractString;
                         time_col::Symbol = :time,
                         default_step_seconds::Real = 3600.0)
    df = CSV.read(csv_path, DataFrame)

    if !(time_col in Symbol.(names(df)))
        nrow(df) > 0 || throw(ArgumentError("Cannot infer time column for an empty schedule CSV."))
        df[!, time_col] = collect(0.0:float(default_step_seconds):(nrow(df) - 1) * float(default_step_seconds))
    end

    _validate_schedule_dataframe(df; time_col = time_col)
    return df
end
