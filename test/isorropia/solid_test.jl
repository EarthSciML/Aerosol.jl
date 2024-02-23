solids = ISORROPIA.generate_solids(t)
@test length(values(solids)) == 19
# Test that the activity of a solid is equal to 1.
@test Symbolics.value(ISORROPIA.activity(solids[:K2SO4])) == 1.0