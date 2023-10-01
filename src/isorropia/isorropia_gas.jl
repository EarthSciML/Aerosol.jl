"""
A gas with the given partial pressure.
"""
struct Gas <: Species
    p
end
"""
From Section 2.2 in Fountoukis and Nenes (2007), the activity of a gas is its partial pressure (in atm).
"""
activity(g::Gas) = g.p
γ(g::Gas) = 1.0
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
        @species $varname($t)=1.0e-6 [bounds=($lbound, $ubound), unit = u"atm", description=$description]
        $varname = ParentScope($varname)
        push!($all_gases, $varname)
        $gasname = $Gas($varname)
    end)
end

# Tests
@test length(all_gases) == 4
@test ModelingToolkit.get_unit(activity(HNO3g)) == u"atm"