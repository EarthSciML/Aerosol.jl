@test length(all_solids) == 19
# Test that the activity of a solid is equal to 1.
@test Symbolics.value(activity(K2SO4s)) == 1.0