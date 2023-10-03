"""
A gas with the given partial pressure.
"""
struct Gas <: Species
    p
end
"""
From Section 2.2 in Fountoukis and Nenes (2007), the activity of a gas is its partial pressure (in atm).
"""
γ(g::Gas) = mol2atm(1.0)
terms(g::Gas) = [g.p], [1]

# Generate the gases.
# Each gas has an associated MTK variable named 
# <name>_g, where <name> is the name of the compound, and
# a Gas struct named <name>g.
all_gases = []
for s ∈ (:HNO3, :HCl, :NH3, :SO4)
    varname = Symbol(s, "_g")
    gasname = Symbol(s, "g")
    description = "Gasous $s"
    eval(quote
        @species $varname($t)=1.0e-6 [unit = u"mol/m_air^3", description=$description]
        $varname = ParentScope($varname)
        push!($all_gases, $varname)
        $gasname = $Gas($varname)
    end)
end

# Tests
@test length(all_gases) == 4
@test ModelingToolkit.get_unit(activity(HNO3g)) == u"atm"

# Test that activity is equal to partial pressum in atm.
@test isequal(ModelingToolkit.subs_constants(activity(HClg)),
            ModelingToolkit.subs_constants(mol2atm(HClg.p)))