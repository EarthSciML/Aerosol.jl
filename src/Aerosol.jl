module Aerosol

using Reexport
using DynamicQuantities
using ModelingToolkit
using ModelingToolkit: t, D
using EarthSciMLBase

@register_unit ppb 1

include("VBS.jl")
include("elemental_carbon.jl")
include("isorropia/isorropia.jl")
@reexport using .ISORROPIA
include("seinfeld_pandis/seinfeld_pandis.jl")
@reexport using .SeinfeldPandis

end
