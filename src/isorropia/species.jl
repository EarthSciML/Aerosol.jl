"""
A species represents a chemical species in the system.

It should have a `terms` method which returns the variable(s) and
stoichiometry coefficient(s) associated with the species, and 
a `vars` method which returns only the variable(s) associated with the species.
"""
abstract type Species end

"Return the combined terms for a vector of species."
function terms(s::AbstractVector) 
   tt = terms.(s)
   vcat([t[1] for t ∈ tt]...), vcat([t[2] for t ∈ tt]...)
end

# Miscellaneous variables and parameters

@parameters T = 293.15 [unit = u"K", description = "Temperature"]
@parameters RH = 0.3 [description = "Relative humidity (expressed on a scale from 0 to 1)"] # unitless
for p ∈ (:T, :RH)
   eval(:($p = ParentScope($p)))
end

#@constants tiny_conc=1e-20 [unit = u"mol/m_air^3", description = "Tiny concentration to avoid division by zero"]

#@constants R = 8.31446261815324 [unit = u"m_air^3*Pa/K/mol", description = "Universal gas constant"]
#@constants PaPerAtm = 101325 [unit = u"Pa/atm", description = "Number of pascals per atmosphere"]
#@constants W_one = 1 [unit = u"kg_water/m_air^3", description = "Unit aerosol water content"]
#mol2atm(p) = p / PaPerAtm * R * T
#atm2mol(p) = p * PaPerAtm / R / T