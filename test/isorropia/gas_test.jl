# Tests
gases = ISORROPIA.generate_gases(t)
@test length(gases) == 3
@test ModelingToolkit.get_unit(ISORROPIA.activity(gases[:HNO3])) == u"atm"

# Test that activity is equal to partial pressum in atm.
@test isequal(ModelingToolkit.subs_constants(ISORROPIA.activity(gases[:HCl])),
    ModelingToolkit.subs_constants(gases[:HCl].p))