"""
A species represents a chemical species in the system.

Chemical species should have an `γ` method which returns the activity coefficient
of the species as specified in Section 2.2 of Fountoukis and Nenes (2007). 
They should also have a `terms` method which returns the variable(s) and
stoichiometry coefficient(s) associated with the species, and also a 
`min_conc` method which returns the minimum concentration among the terms that 
make up the species.
"""
abstract type Species end
activity(s) = reduce(*, [m^ν for (m, ν) ∈ zip(terms(s)...)]) * γ(s)
γ(s::Species) = error("activity coefficient γ not defined for $(typeof(s))")
terms(s::Species) = error("terms not defined for $(typeof(s))")

""" 
The activity coefficient of multiple species is the product of their Activity
coefficients as shown in Table 2 of Fountoukis and Nenes (2007).
"""
γ(s::AbstractVector) = reduce(*, γ.(s))
function terms(s::AbstractVector) 
    tt = terms.(s)
    vcat([t[1] for t ∈ tt]...), vcat([t[2] for t ∈ tt]...)
end
min_conc(s::AbstractArray) = reduce(min, min_conc.(s))

# Miscellaneous variables and parameters
@variables I(t) = 1.0e-4 [unit = u"mol/kg_water", description = "Ionic strength"]
@variables W(t) = 1.0e-8 [unit = u"kg_water/m_air^3", description = "Aerosol water content"]
for v ∈ (:I, :W)
    eval(:($v = ParentScope($v))) # Keep these as global MTK variables.
end

@parameters T = 293.15 [unit = u"K", description = "Temperature"]
@parameters RH = 0.3 [description = "Relative humidity (expressed on a scale from 0 to 1)"] # unitless
for p ∈ (:T, :RH)
   eval(:($p = ParentScope($p)))
end

@constants tiny_conc=1e-20 [unit = u"mol/m_air^3", description = "Tiny concentration to avoid division by zero"]

@constants R = 8.31446261815324 [unit = u"m_air^3*Pa/K/mol", description = "Universal gas constant"]
@constants PaPerAtm = 101325 [unit = u"Pa/atm", description = "Number of pascals per atmosphere"]
@constants W_one = 1 [unit = u"kg_water/m_air^3", description = "Unit aerosol water content"]
mol2atm(p) = p / PaPerAtm * R * T