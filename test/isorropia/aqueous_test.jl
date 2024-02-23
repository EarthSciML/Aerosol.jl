@variables t [description = "Time" unit = u"s"]

ions = ISORROPIA.generate_ions(t)
salts = ISORROPIA.generate_salts(ions)
@test length(ions) == 14
@test length(salts) == 23

active_salts = collect(values(salts))
@test issetequal(ISORROPIA.same_cation(salts[:KCl], active_salts, salts),
    [salts[:K2SO4], salts[:KNO3], salts[:KCl], salts[:KHSO4]])
@test issetequal(ISORROPIA.same_anion(salts[:KCl], active_salts, salts),
    [salts[:CaCl2], salts[:KCl], salts[:MgCl2], salts[:NaCl], salts[:NH4Cl], salts[:HCl]])

W = ISORROPIA.Water(t, active_salts)
I = ISORROPIA.IonicStrength(t, values(ions), W.W)

@test ModelingToolkit.get_unit(ISORROPIA.logγ₁₂T⁰(salts[:NaCl], active_salts, salts, I.I, W.W)) == u"mol^0.5/kg_water^0.5"
@test ModelingToolkit.get_unit(ISORROPIA.F₁(salts[:NaCl], active_salts, salts, I.I, W.W)) == u"mol^0.5/kg_water^0.5"
@test ModelingToolkit.get_unit(ISORROPIA.F₂(salts[:NaCl], active_salts, salts, I.I, W.W)) == u"mol^0.5/kg_water^0.5"
@test ModelingToolkit.get_unit(ISORROPIA.X(salts[:NaCl], I.I, W.W)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(ISORROPIA.Y(salts[:NaCl], I.I, W.W)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(ISORROPIA.Γ⁰(salts[:NaCl].q, I.I)) == u"mol/kg_water"
@test ModelingToolkit.get_unit(ISORROPIA.Γstar(1, I.I)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(ISORROPIA.C(1, I.I)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(ISORROPIA.logγ⁰₁₂(salts[:NaCl], I.I)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(ISORROPIA.logγ⁰₁₂(salts[:NaCl], I.I)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(ISORROPIA.logγ₁₂(salts[:NaCl], active_salts, salts, I.I, W.W)) isa Unitful.FreeUnits{(),NoDims,nothing}

function sub(expr, u=nothing)
    if !isnothing(u)
        expr = substitute(expr, u)
    end
    expr = ModelingToolkit.subs_constants(expr)
    defaults = [I.I => 2.5, ions[:Cl].m => 1, ions[:H].m => 2, ISORROPIA.T => 289,
        ions[:K].m => 0.5, ions[:Mg].m => 1.2, ions[:NH4].m => 2.5, 
        ions[:NO3].m => 3.1, ions[:Na].m => 0.2, ions[:Ca].m => 1.2,
        ions[:SO4].m => 2.0, ions[:HSO4].m => 0.8, W.W => 0.8]
    substitute(expr, defaults)
end

# Activity coefficients should often be ≤ 1 for ions in aqueous solution
@test sub(ISORROPIA.logγ⁰₁₂(salts[:KCl], I.I)) ≈ -0.16013845145909214

@test sub(ISORROPIA.logγ₁₂T⁰(salts[:KCl], active_salts, salts, I.I, W.W)) ≈ 0.8322034507224931

@test sub(ISORROPIA.logγ₁₂(salts[:KCl], active_salts, salts, I.I, W.W)) ≈ 0.8859124609238154

# Activity coefficients should decrease with increasing ionic strength.
@test sub(ISORROPIA.logγ₁₂(salts[:KCl], active_salts, salts, I.I, W.W), [I.I => 5]) ≈ 0.853546306554723

# Test activity coefficients for all salts.
saltnames = [:CaNO32, :CaCl2, :K2SO4, :KNO3, :KCl, :MgSO4,
    :MgNO32, :MgCl2, :NaCl, :Na2SO4, :NaNO3, :NH42SO4, :NH4NO3,
    :NH4Cl, :H2SO4, :HHSO4, :HNO3, :HCl]
want_vals = [0.9899740297791452, 2.219122534658643, -1.2832826882588682, -0.03594891773580812, 0.8859124609238154,
    0.6655782776789243, 1.6790696497304067, 2.9082181546099055, 1.288358433201964, -0.7466880585546696, 0.3664970545423407,
    -1.0102747960911531, 0.16880700138997853, 1.0906683800496015, 0.1681542434957075, 0.7541249743955618, 1.0526287810801234, 1.9744901597397468]
for (i, salt) in enumerate(saltnames)
    v = sub(ISORROPIA.logγ₁₂(salts[salt], active_salts, salts, I.I, W.W))
    @test v ≈ want_vals[i]
end

# Units in last column in Table 2.

@test ModelingToolkit.get_unit(ISORROPIA.activity(salts[:NaCl], active_salts, salts, I.I, W.W)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(ISORROPIA.activity(salts[:CaNO32], active_salts, salts, I.I, W.W)) == u"mol^3/kg_water^3"

@test sub(ISORROPIA.activity(salts[:NaCl], active_salts, salts, I.I, W.W)) ≈ 4.110587887491537

# Special cases

# Units in last column in Table 2.
unitactf(x) = ModelingToolkit.get_unit(ISORROPIA.activity(salts[x], active_salts, salts, I.I, W.W)) 
@test unitactf(:CaSO4) == u"mol^2/kg_water^2"
@test unitactf(:KHSO4) == u"mol^2/kg_water^2"
@test unitactf(:NH4HSO4) == u"mol^2/kg_water^2"
@test unitactf(:NaHSO4) == u"mol^2/kg_water^2"
@test unitactf(:NH43HSO42) == u"mol^5/kg_water^5"

actf(x) = sub(ISORROPIA.activity(salts[x], active_salts, salts, I.I, W.W)) 
want_γ = Any[1.0e-20, 0.8460080854008434, 0.9372095310924331, 1.0345811139028118, 0.5384101987210793]
want_activity = Any[3.7499999999999995e-40, 0.4473310503522503, 2.744880328657807, 0.2675895203110956, 1.3807544582351838]
for (i, s) ∈ enumerate([:CaSO4, :KHSO4, :NH4HSO4, :NaHSO4, :NH43HSO42])
    @test sub(ISORROPIA.γ₁₂(salts[s], active_salts, salts, I.I, W.W)) ≈ want_γ[i]
    @test actf(s) ≈ want_activity[i]
end

@test nameof(salts[:CaCl2]) == "CaCl2_aqs"