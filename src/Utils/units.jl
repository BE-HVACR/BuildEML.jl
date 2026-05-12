using Unitful


to_SI(x, u) = x isa Unitful.Quantity ? Float64(Unitful.ustrip(Unitful.uconvert(u, x))) : Float64(x)

# Common

to_SI_T(x) = x isa Unitful.Quantity ? Float64(Unitful.ustrip(Unitful.uconvert(u"K", x))) : Float64(x)

to_SI_p(x) = x isa Unitful.Quantity ? Float64(Unitful.ustrip(Unitful.uconvert(u"Pa", x))) : Float64(x)

to_SI_massflow(x) = x isa Unitful.Quantity ? Float64(Unitful.ustrip(Unitful.uconvert(u"kg*s^-1", x))) : Float64(x)

to_SI_volflow(x) = x isa Unitful.Quantity ? Float64(Unitful.ustrip(Unitful.uconvert(u"m^3*s^-1", x))) : Float64(x)


to_SI_energy(x) = x isa Unitful.Quantity ? Float64(Unitful.ustrip(Unitful.uconvert(u"J", x))) : Float64(x)

to_SI_power(x) = x isa Unitful.Quantity ? Float64(Unitful.ustrip(Unitful.uconvert(u"W", x))) : Float64(x)

unit_K2C(TK) = TK - 273.15
unit_C2K(TC) = TC + 273.15


function ensure_SI(; T=nothing, p=nothing, ṁ=nothing, V̇=nothing, E=nothing, P=nothing)
    return (; 
        T = isnothing(T) ? nothing : to_SI_T(T),
        p = isnothing(p) ? nothing : to_SI_p(p),
        ṁ = isnothing(ṁ) ? nothing : to_SI_massflow(ṁ),
        V̇ = isnothing(V̇) ? nothing : to_SI_volflow(V̇),
        E = isnothing(E) ? nothing : to_SI_energy(E),
        P = isnothing(P) ? nothing : to_SI_power(P),
    )
end


