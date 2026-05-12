# Weather file lookup for bundled EPW files.
#
# Public API:
#   weather_file_path(city; file_type="epw")            → path by city name or alias (case-insensitive)
#   typical_climate_weather_path(zone; file_type="epw") → path by ASHRAE climate zone (e.g. "5A")
#
# City data currently covers ASHRAE 169-2025 climate zones 1A–8 plus Miami and Chicago.
# EPW files live in the weatherfile/ subdirectory next to this file.

const WEATHERFILE_DIR = normpath(joinpath(@__DIR__, "weatherfile"))

const WEATHER_CITY_DATA = [
    (
        canonical = "Honolulu",
        climate_zone = "1A",
        aliases = ["Honolulu"],
        files = Dict("epw" => "USA_HI_CGB.Honolulu.Oahu.994007_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Yuma",
        climate_zone = "1B",
        aliases = ["Yuma", "Yuma Airport", "Yuma MCAS", "Phoenix", "Phoenix Sky Harbor", "Sky Harbor"],
        files = Dict("epw" => "USA_AZ_Yuma-MCAS.Yuma-Yuma.Intl.AP.722800_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Houston",
        climate_zone = "2A",
        aliases = ["Houston"],
        files = Dict("epw" => "USA_TX_Houston-Bush.Intercontinental.AP.722430_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Tucson",
        climate_zone = "2B",
        aliases = ["Tucson"],
        files = Dict("epw" => "USA_AZ_Tucson.722740_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Atlanta",
        climate_zone = "3A",
        aliases = ["Atlanta"],
        files = Dict("epw" => "USA_GA_Atlanta-Hartsfield-Jackson.Intl.AP.722190_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "El Paso",
        climate_zone = "3B",
        aliases = ["El Paso", "ElPaso"],
        files = Dict("epw" => "USA_TX_El.Paso.Intl.AP.722700_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "San Francisco",
        climate_zone = "3C",
        aliases = ["San Francisco", "SanFrancisco", "SFO"],
        files = Dict("epw" => "USA_CA_San.Francisco.Intl.AP.724940_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "New York",
        climate_zone = "4A",
        aliases = ["New York", "NewYork", "JFK"],
        files = Dict("epw" => "USA_NY_New.York-Kennedy.Intl.AP.744860_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Albuquerque",
        climate_zone = "4B",
        aliases = ["Albuquerque"],
        files = Dict("epw" => "USA_NM_Albuquerque.Intl.Sunport.723650_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Seattle",
        climate_zone = "4C",
        aliases = ["Seattle"],
        files = Dict("epw" => "USA_WA_Seattle-Tacoma.Intl.AP.727930_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Buffalo",
        climate_zone = "5A",
        aliases = ["Buffalo"],
        files = Dict("epw" => "USA_NY_Buffalo.Niagara.Intl.AP.725280_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Kit Carson",
        climate_zone = "5B",
        aliases = ["Kit Carson", "KitCarson", "Burlington", "Burlington Kit Carson", "Denver"],
        files = Dict("epw" => "USA_CO_Burlington-Kit.Carson.County.AP.724689_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Port Angeles",
        climate_zone = "5C",
        aliases = ["Port Angeles", "PortAngeles"],
        files = Dict("epw" => "USA_WA_Port.Angeles-Fairchild.Intl.AP.727885_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Rochester",
        climate_zone = "6A",
        aliases = ["Rochester"],
        files = Dict("epw" => "USA_MN_Rochester.Intl.AP.726440_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Great Falls",
        climate_zone = "6B",
        aliases = ["Great Falls", "GreatFalls"],
        files = Dict("epw" => "USA_MT_Great.Falls.727760_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "International Falls",
        climate_zone = "7",
        aliases = ["International Falls", "InternationalFalls"],
        files = Dict("epw" => "USA_MN_International.Falls-Falls.Intl.AP.727470_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Fairbanks",
        climate_zone = "8",
        aliases = ["Fairbanks"],
        files = Dict("epw" => "USA_AK_Fairbanks.Intl.AP.702610_TMYx.2011-2025.epw"),
    ),
    (
        canonical = "Miami",
        climate_zone = nothing,
        aliases = ["Miami"],
        files = Dict("epw" => "USA_FL_Miami.Intl.AP.722020_TMY3.epw"),
    ),
    (
        canonical = "Chicago",
        climate_zone = nothing,
        aliases = ["Chicago", "OHare", "O'Hare"],
        files = Dict("epw" => "USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.epw",
                     "csv" => "USA_IL_Chicago-OHare.Intl.AP.725300_TMY3EPW.csv",
                     "ddy" => "USA_IL_Chicago-OHare.Intl.AP.725300_TMY3EPW.ddy",
                     "stat" => "USA_IL_Chicago-OHare.Intl.AP.725300_TMY3EPW.stat"),
    ),
]

_normalize_lookup_key(x)   = replace(lowercase(strip(String(x))), r"[^a-z0-9]+" => "")
_normalize_file_type(x)    = replace(lowercase(strip(String(x))), "." => "")
_normalize_climate_zone(x) = replace(uppercase(strip(String(x))), r"[^A-Z0-9]+" => "")

const WEATHER_CITY_FILES = Dict(item.canonical => item.files for item in WEATHER_CITY_DATA)

const _WEATHER_CITY_ALIASES = let aliases = Dict{String,String}()
    for item in WEATHER_CITY_DATA
        aliases[_normalize_lookup_key(item.canonical)] = item.canonical
        for alias in item.aliases
            aliases[_normalize_lookup_key(alias)] = item.canonical
        end
    end
    aliases
end

const ASHRAE_TYPICAL_CLIMATE_CITY = Dict(
    item.climate_zone => item.canonical for item in WEATHER_CITY_DATA if !isnothing(item.climate_zone)
)

const _ASHRAE_TYPICAL_CLIMATE_CITY_ALIASES = Dict(
    _normalize_climate_zone(zone) => city for (zone, city) in ASHRAE_TYPICAL_CLIMATE_CITY
)

function _canonical_weather_city(city)
    key = _normalize_lookup_key(city)
    canonical = get(_WEATHER_CITY_ALIASES, key, nothing)
    canonical === nothing && throw(ArgumentError("Unknown weather city $(repr(city)). Available cities: $(join(sort!(collect(keys(WEATHER_CITY_FILES))), ", "))"))
    return canonical
end

function _typical_climate_city(zone)
    key = _normalize_climate_zone(zone)
    city = get(_ASHRAE_TYPICAL_CLIMATE_CITY_ALIASES, key, nothing)
    city === nothing && throw(ArgumentError("Unknown ASHRAE climate zone $(repr(zone)). Available zones: $(join(sort!(collect(keys(ASHRAE_TYPICAL_CLIMATE_CITY))), ", "))"))
    return city
end

function weather_file_path(city; file_type::Union{AbstractString,Symbol} = "epw")
    canonical = _canonical_weather_city(city)
    ext = _normalize_file_type(file_type)
    files = WEATHER_CITY_FILES[canonical]
    filename = get(files, ext, nothing)
    if filename === nothing
        available = join(sort!(collect(keys(files))), ", ")
        throw(ArgumentError("Weather file type $(repr(file_type)) is not available for $(canonical). Available types: $(available)"))
    end
    path = normpath(joinpath(WEATHERFILE_DIR, filename))
    isfile(path) || throw(ArgumentError("Mapped weather file does not exist on disk: $(path)"))
    return path
end

function typical_climate_weather_path(zone; file_type::Union{AbstractString,Symbol} = "epw")
    return weather_file_path(_typical_climate_city(zone); file_type = file_type)
end
