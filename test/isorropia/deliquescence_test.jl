@variables RH [description = "Relative Humidity"]

ions = ISORROPIA.generate_ions(t)
salts = ISORROPIA.generate_salts(ions)
active_salts = collect(values(salts))

@test ISORROPIA.drh(salts[:KNO3]) == salts[:KNO3].drh

@test ModelingToolkit.substitute(ModelingToolkit.subs_constants(ISORROPIA.drh(salts[:CaNO32])), ISORROPIA.T => 298.15) == salts[:CaNO32].drh

@test ModelingToolkit.substitute(ModelingToolkit.subs_constants(ISORROPIA.drh(salts[:CaNO32])), ISORROPIA.T => 320) ≈ 0.5513060522349494

@test ModelingToolkit.get_unit(ISORROPIA.drh(salts[:CaNO32])) isa Unitful.FreeUnits{(),NoDims,nothing}

# TODO(CT): Our solution MDRH selection doesn't work in most cases, 
# because our method of checking which ions are present doesn't 
# yield unique matches. We need to find a better 
# way to do this, perhaps based on the ratios in Fountoukis and 
# Nenes (2007) Table 3. 
@testset "solution_mdrh_recurrent" begin
    for i ∈ eachindex(ISORROPIA.mdrhs)
        u = Dict()
        for ion ∈ values(ions)
            u[ion.m] = 1.e-20
            u[ion.m] = 1.e-20
        end
        for s ∈ ISORROPIA.mdrhs[i][1]
            u[salts[s].cation.m] = 1.e-9
            u[salts[s].anion.m] = 1.e-9
        end
        @testset "$i" begin
            x = ModelingToolkit.substitute(ModelingToolkit.subs_constants(ISORROPIA.solution_mdrh_recurrent(1, active_salts, salts, ions)), u)
            if i ∈ [3, 5, 9, 10, 11, 12, 13, 14]
                @test_broken x == ISORROPIA.mdrhs[i][2]
            else 
                @test x == ISORROPIA.mdrhs[i][2]
            end 
        end
    end
end

@testset "f_drhs" begin
    del = ISORROPIA.deliquescence(t, RH, active_salts, salts, ions)[1]
    @unpack DRH_NH43HSO42_aqs, f_NH43HSO42_aqs, MDRH = del
    index = [isequal(eq.lhs, f_NH43HSO42_aqs) for eq in equations(del)]

    subrh(rh) = begin
        subbed = substitute(del, [RH => rh, defaults...])
        only(equations(subbed)[index]).rhs
    end

    defaults = [del.metastable => 0.0, DRH_NH43HSO42_aqs => 0.7, MDRH => 0.4]
    @test subrh(0.71) == 0.0
    @test subrh(0.39) == 1.0
    @test subrh(0.55) ≈ 0.5
end