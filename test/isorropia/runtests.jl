using Unitful, ModelingToolkit, Catalyst
using Test

@variables t [unit = u"s", description = "Time"]

include(joinpath(@__DIR__, "../../src/isorropia/units.jl"))
include(joinpath(@__DIR__, "../../src/isorropia/species.jl"))
include(joinpath(@__DIR__, "../../src/isorropia/aqueous.jl"))

@testset "aqueous" begin include("aqueous_test.jl") end

include(joinpath(@__DIR__, "../../src/isorropia/deliquescence.jl"))

@testset "deliquescence" begin include("deliquescence_test.jl") end

include(joinpath(@__DIR__, "../../src/isorropia/gas.jl"))

@testset "gas" begin include("gas_test.jl") end

include(joinpath(@__DIR__, "../../src/isorropia/solid.jl"))

@testset "solid" begin include("solid_test.jl") end

include(joinpath(@__DIR__, "../../src/isorropia/water.jl"))

@testset "water" begin include("water_test.jl") end

include(joinpath(@__DIR__, "../../src/isorropia/reactions.jl"))

@testset "reactions" begin include("reactions_test.jl") end

@safetestset "isorropia" begin include("isorropia_test.jl") end