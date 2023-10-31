"""
A solid with the given concentration.
"""
struct Solid <: Species
    m
end

#@constants unit_conc = 1.0 [unit = u"mol/m_air^3", description = "Unit concentration"]

"""
From Section 2.2 in Fountoukis and Nenes (2007), the activity of a solid is 1.
"""
activity(s::Solid) = 1.0

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
