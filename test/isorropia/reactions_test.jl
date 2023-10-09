@test ModelingToolkit.get_unit(mol2atm(SO4_g)) == u"atm"

# Units from Table 2, last column
@test ModelingToolkit.get_unit(rxn5.sys.K_eq) == u"mol^2/kg_water^2"

# Test double activity
@test ModelingToolkit.get_unit(activity([NH3g, HNO3g])) == u"atm^2"
