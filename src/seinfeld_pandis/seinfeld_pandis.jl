"""
Aerosol thermodynamics equations from Seinfeld & Pandis (2006) Chapter 10.

Reference: Seinfeld, J.H. and Pandis, S.N. (2006) "Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change", 2nd Edition, John Wiley & Sons, Chapter 10.
"""
module SeinfeldPandis

using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities

export KelvinEffect, DRHTemperature, ZSRWaterContent, NH4NO3Equilibrium
export kelvin_saturation_ratio, drh_temperature, nh4no3_drh, zsr_water_content
export nh4no3_Kp, nh4no3_K_AN, ionic_strength_fraction, get_equilibrium_constant

include("kelvin_effect.jl")
include("drh.jl")
include("zsr.jl")
include("nh4no3_equilibrium.jl")

end # module
