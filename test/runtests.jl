using Aerosol
using Test
using SafeTestsets

@testset "Aerosol.jl" begin
    @safetestset "isorropia" begin include(joinpath(@__DIR__, "isorropia", "runtests.jl")) end
end
