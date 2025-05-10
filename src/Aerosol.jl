module Aerosol

using Reexport
using DynamicQuantities
using ModelingToolkit
using ModelingToolkit: t, D

@register_unit ppb 1

#include("isorropia/isorropia.jl")
#@reexport using .ISORROPIA
#@reexport using .ISORROPIA.MyUnits
include("elemental_carbon.jl")

end
