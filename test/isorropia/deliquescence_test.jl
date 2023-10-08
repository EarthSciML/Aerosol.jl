@test drh(KNO3_aqs) == KNO3_aqs.drh

@test ModelingToolkit.substitute(ModelingToolkit.subs_constants(drh(CaNO32_aqs)), T => 298.15) == CaNO32_aqs.drh

@test ModelingToolkit.substitute(ModelingToolkit.subs_constants(drh(CaNO32_aqs)), T => 320) ≈ 0.5513060522349494

@test ModelingToolkit.get_unit(drh(CaNO32_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}

# TODO(CT): Our solution MDRH selection doesn't work in most cases, 
# because our method of checking which ions are present doesn't 
# yield unique matches. We need to find a better 
# way to do this, perhaps based on the ratios in Fountoukis and 
# Nenes (2007) Table 3. 
@testset "solution_mdrh_recurrent" begin
    for i ∈ eachindex(mdrhs)
        u = Dict()
        for ion ∈ all_ions
            u[ion] = 1.e-20
            u[ion] = 1.e-20
        end
        for s ∈ mdrhs[i][1]
            u[s.cation.m] = 1.e-10
            u[s.anion.m] = 1.e-10
        end
        @testset "$i" begin
            x = ModelingToolkit.substitute(ModelingToolkit.subs_constants(solution_mdrh_recurrent(1)), u)
            if i ∈ [3, 5, 9, 10, 11, 12, 13, 14]
                @test_broken x == mdrhs[i][2]
            else 
                @test x == mdrhs[i][2]
            end 
        end
    end
end


@testset "f_drhs" begin
    defaults = [metastable => 0.0, DRH_NH43HSO42 => 0.7, MDRH => 0.4]
    @test substitute(drh_eqs[end].rhs, [RH => 0.71, defaults...]) == 0.0
    @test substitute(drh_eqs[end].rhs, [RH => 0.39, defaults...]) == 1.0
    @test substitute(drh_eqs[end].rhs, [RH => 0.55, defaults...]) ≈ 0.5
end