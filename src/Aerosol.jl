module Aerosol

using Reexport
using DynamicQuantities
using ModelingToolkit
using ModelingToolkit: t, D
using EarthSciMLBase

@register_unit ppb 1

include("VBS.jl")
include("elemental_carbon.jl")
include("single_particle_dynamics.jl")
include("isorropia/isorropia.jl")
@reexport using .ISORROPIA

# Aqueous chemistry module (Seinfeld & Pandis Chapter 7)
include("henrys_law.jl")
include("aqueous_equilibria.jl")
include("sulfate_formation.jl")
include("cloud_chemistry.jl")

# Export aqueous chemistry components
export HenrysLaw, HenrysLawTemperature, EffectiveHenrysLaw
export WaterEquilibrium, CO2Equilibria, SO2Equilibria, NH3Equilibria
export HNO3Equilibria, H2O2Equilibria, O3Equilibria, AqueousEquilibria
export SulfateFormationO3, SulfateFormationH2O2
export SulfateFormationFe, SulfateFormationMn, SulfateFormationFeMn
export SulfateFormation
export CloudChemistry, CloudChemistryFixedpH, CloudChemistryODE

# Export aqueous chemistry utility functions
export effective_henrys_constant, aqueous_fraction, distribution_factor
export henrys_constant_at_T
export rate_to_ppb_hr, rate_to_percent_hr, so2_lifetime

# Export aqueous chemistry constants
export ATM_TO_PA, R_GAS_PA, R_GAS_J
export HENRY_CONSTANTS_298, DELTA_H_DISSOLUTION
export K_W_298, K_C1_298, K_C2_298, K_S1_298, K_S2_298, K_A1_298
export K0_O3, K1_O3, K2_O3, K_H2O2, K_MN_SIV, K_FE_SIV, K_FEMN_SIV

end
