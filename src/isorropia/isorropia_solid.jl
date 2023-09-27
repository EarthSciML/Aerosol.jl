"""
A solid with the given concentration.
"""
struct Solid <: Species
    m
end
"""
From Section 2.2 in Fountoukis and Nenes (2007), the activity of a solid is 1.
"""
activity(s::Solid) = 1.0

# Generate the solid compounds.
# Each solid compound has an associated MTK variable named 
# <name>_s, where <name> is the name of the compound, and
# a Solid struct named <name>s.
all_solids = []
for s = (:CaNO32, :CaCl2, :CaSO4, :KHSO4, :K2SO4, :KNO3, :KCl, 
    :MgSO4, :MgNO32, :MgCl2, :NaCl, :NaNO3, :Na2SO4, :NaHSO4, :NH4Cl, 
    :NH4NO3, :NH42SO4, :NH4HSO4, :NH43HSO42)
    varname = Symbol(s, "_s")
    solidname = Symbol(s, "s")
    description = "Solid $s"
    eval(quote
        @variables $varname=0 [unit = u"mol/m_air^3", description=$description]
        push!($all_solids, $varname)
        $solidname = $Solid($varname)
    end)
end

# Tests
@test length(all_solids) == 19
@test activity(K2SO4s) === 1.0