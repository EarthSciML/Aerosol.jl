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

"""
Generate the gases.
Each gas has an associated MTK variable named <name>_g, where <name> is the name of the compound.
All SO4 immediately goes to aerosol phase as per Section 3.3 (item 1) of Fountoukis and Nenes (2007), so we don't include it here.
"""
function generate_gases(t)
    gases = Dict()
    for (s, v) ∈ [:HNO3 => 1e-8, :HCl => 1e-8, :NH3 => 1e-8]
        varname = Symbol(s, "_g")
        description = "Gasous $s"
        var = only(@species $varname(t) = v [unit = u"Constants.atm", description = description])
        var = ParentScope(var)
        gases[s] = Gas(var)
    end
    gases
end
