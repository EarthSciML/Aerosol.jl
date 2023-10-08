# Tests
@test length(all_gases) == 4
@test ModelingToolkit.get_unit(activity(HNO3g)) == u"atm"

# Test that activity is equal to partial pressum in atm.
@test isequal(ModelingToolkit.subs_constants(activity(HClg)),
    ModelingToolkit.subs_constants(mol2atm(HClg.p)))