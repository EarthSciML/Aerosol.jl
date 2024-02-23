using Unitful, ModelingToolkit, Catalyst
using Test
using EarthSciMLBase
using DifferentialEquations
using Plots

include(joinpath(@__DIR__, "../../src/isorropia/isorropia.jl"))
using .ISORROPIA
using Main.ISORROPIA.MyUnits

@variables t [unit = u"s", description = "Time"]

@testset "aqueous" begin include("aqueous_test.jl") end

@testset "deliquescence" begin include("deliquescence_test.jl") end

@testset "gas" begin include("gas_test.jl") end

@testset "solid" begin include("solid_test.jl") end

@testset "water" begin include("water_test.jl") end

@testset "reactions" begin include("reactions_test.jl") end

@safetestset "isorropia" begin include("isorropia_test.jl") end