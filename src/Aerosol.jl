module Aerosol

using Reexport
using DynamicQuantities
using ModelingToolkit
using ModelingToolkit: t, D
using EarthSciMLBase

@register_unit ppb 1

include("VBS.jl")
include("elemental_carbon.jl")
include("mass_transfer.jl")
include("timescales.jl")
include("aqueous_transport.jl")
include("isorropia/isorropia.jl")
@reexport using .ISORROPIA

end
