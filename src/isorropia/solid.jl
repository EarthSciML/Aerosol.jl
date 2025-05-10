"""
A solid with the given concentration.
"""
struct Solid <: Species
    m
end

"""
From Section 2.2 in Fountoukis and Nenes (2007), the activity of a solid is 1. 
However, this cannot be strictly true because if the activity were always equal to one the solid would continue to be consumed 
even after its concentration had reached zero, resulting in negative concentrations which are physically impossible.
"""
activity(s::Solid) = 1

vars(s::Solid) = [s.m]
terms(s::Solid) = [s.m], [1]

function Base.nameof(s::Solid)
    string(Symbolics.tosymbol(s.m, escape=false))
end

"""
Generate the solid compounds.
Each solid compound has an associated MTK variable named 
<name>_s, where <name> is the name of the compound, and
a Solid struct named <name>s.
"""
function generate_solids(t)
    solids = Dict()
    for s = (:CaNO32, :CaCl2, :CaSO4, :KHSO4, :K2SO4, :KNO3, :KCl,
        :MgSO4, :MgNO32, :MgCl2, :NaCl, :NaNO3, :Na2SO4, :NaHSO4, :NH4Cl,
        :NH4NO3, :NH42SO4, :NH4HSO4, :NH43HSO42)
        varname = Symbol(s, "_s")
        description = "Solid $s"
        v = only(@variables $varname(t) = 1e-8 [unit = u"mol/m_air^3", description = description])
        v = ParentScope(v)
        solids[s] = Solid(v)
    end
    solids
end
