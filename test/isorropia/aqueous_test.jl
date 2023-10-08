@test length(all_ions) == 14
@test length(all_salts) == 23

@test same_cation(KCl_aqs) == [K2SO4_aqs, KNO3_aqs, KCl_aqs, KHSO4_aqs]
@test same_anion(KCl_aqs) == [CaCl2_aqs, KCl_aqs, MgCl2_aqs, NaCl_aqs, NH4Cl_aqs, HCl_aqs]

@test ModelingToolkit.get_unit(logγ₁₂T⁰(NaCl_aqs)) == u"mol^0.5/kg_water^0.5"
@test ModelingToolkit.get_unit(F₁(NaCl_aqs)) == u"mol^0.5/kg_water^0.5"
@test ModelingToolkit.get_unit(F₂(NaCl_aqs)) == u"mol^0.5/kg_water^0.5"
@test ModelingToolkit.get_unit(X(NaCl_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(Y(NaCl_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(Γ⁰(NaCl_aqs.q)) == u"mol/kg_water"
@test ModelingToolkit.get_unit(Γstar(1)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(C(1)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(logγ⁰₁₂(NaCl_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(logγ⁰₁₂(NaCl_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(logγ₁₂(NaCl_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}

function sub(expr, u=nothing)
    if !isnothing(u)
        expr = substitute(expr, u)
    end
    expr = ModelingToolkit.subs_constants(expr)
    defaults = [I => 2.5, Cl_aq => 1, H_aq => 2, T => 289,
        K_aq => 0.5, Mg_aq => 1.2, NH4_aq => 2.5, NO3_aq => 3.1, Na_aq => 0.2, Ca_aq => 1.2,
        SO4_aq => 2.0, HSO4_aq => 0.8, W => 0.8]
    substitute(expr, defaults)
end

# Activity coefficients should be ≤ 1 for ions in aqueous solution
@test sub(logγ⁰₁₂(KCl_aqs)) ≈ -0.16013845145909214

@test sub(logγ₁₂T⁰(KCl_aqs)) ≈ 0.8322034507224931

# Activity coefficients should decrease with increasing temperature.
@test sub(logγ₁₂(KCl_aqs)) ≈ 0.8859124609238154

@test sub(logγ₁₂(KCl_aqs)) ≈ 0.8859124609238154

# Activity coefficients should decrease with increasing ionic strength.
@test sub(logγ₁₂(KCl_aqs), [I => 5]) ≈ 0.853546306554723

# Test activity coefficients for all salts.
want_vals = [0.9899740297791452, 2.219122534658643, -1.2832826882588682, -0.03594891773580812, 0.8859124609238154,
    0.6655782776789243, 1.6790696497304067, 2.9082181546099055, 1.288358433201964, -0.7466880585546696, 0.3664970545423407,
    -1.0102747960911531, 0.16880700138997853, 1.0906683800496015, 0.1681542434957075, 0.7541249743955618, 1.0526287810801234, 1.9744901597397468]
for (i, salt) in enumerate(all_salts)
    if typeof(salt) <: SpecialSalt
        continue
    end
    v = sub(logγ₁₂(salt))
    @test v ≈ want_vals[i]
end

# Units in last column in Table 2.
@test ModelingToolkit.get_unit(activity(NaCl_aqs)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(activity(CaNO32_aqs)) == u"mol^3/kg_water^3"

@test sub(activity(NaCl_aqs)) ≈ 4.110587887491537

# Special cases

# Units in last column in Table 2.
@test ModelingToolkit.get_unit(activity(CaSO4_aqs)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(activity(KHSO4_aqs)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(activity(NH4HSO4_aqs)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(activity(NaHSO4_aqs)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(activity(NH43HSO42_aqs)) == u"mol^5/kg_water^5"

want_γ = Any[1.0e-20, 0.8460080854008434, 0.9372095310924331, 1.0345811139028118, 0.5384101987210793]
want_activity = Any[3.7499999999999995e-40, 0.4473310503522503, 2.744880328657807, 0.2675895203110956, 1.3807544582351838]
for (i, s) ∈ enumerate([CaSO4_aqs, KHSO4_aqs, NH4HSO4_aqs, NaHSO4_aqs, NH43HSO42_aqs])
    @test sub(γ₁₂(s)) ≈ want_γ[i]
    @test sub(activity(s)) ≈ want_activity[i]
end

@test nameof(CaCl2_aqs) == "CaCl2"