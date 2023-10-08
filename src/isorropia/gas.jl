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
min_conc(g::Gas) = g.p

# Generate the gases.
# Each gas has an associated MTK variable named 
# <name>_g, where <name> is the name of the compound, and
# a Gas struct named <name>g.
all_gases = []
for (s, v) ∈ [:HNO3 => 1e-20, :HCl => 1e-20, :NH3 => 1e-20, :SO4 => 1e-20]
    varname = Symbol(s, "_g")
    gasname = Symbol(s, "g")
    description = "Gasous $s"
    eval(quote
        @species $varname($t) = $v [unit = u"mol/m_air^3", description = $description]
        $varname = ParentScope($varname)
        push!($all_gases, $varname)
        $gasname = $Gas($varname)
    end)
end
