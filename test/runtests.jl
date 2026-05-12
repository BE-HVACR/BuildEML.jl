using BuildEML
using Test

@testset "BuildEML.jl" begin
    @test isdefined(BuildEML, :Utils)
    @test isdefined(BuildEML, :Media)
    @test isdefined(BuildEML, :FundamentalComponents)
    @test isdefined(BuildEML, :HVACComponents)
    @test isdefined(BuildEML, :Disturbances)
    @test isdefined(BuildEML, :BuildingComponents)
end
