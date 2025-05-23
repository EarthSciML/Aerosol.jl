using Aerosol
using Test
using SafeTestsets

@testset "Aerosol.jl" begin
    @safetestset "elemental carbon" begin
        include("elemental_carbon_test.jl")
    end
    @safetestset "isorropia" begin
        include(joinpath(@__DIR__, "isorropia", "runtests.jl"))
    end
    @safetestset "VBS" begin
        include("VBS_test.jl")
    end
end
