using Aerosol
import Aerosol.ISORROPIA: Aqueous, Ion, Salt
using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using NonlinearSolve
using OrdinaryDiffEq
using Test

@testitem "Pivot Salts" begin
    #=
    This test checks for the set of salts that are linearly independent of each other.
    =#
    using LinearAlgebra

    salts = [0 0 0 1 0 0 0 2 0 0 0
             0 0 0 1 0 0 2 0 0 0 0
             0 0 0 1 0 0 0 0 1 0 0
             0 0 0 0 1 0 0 0 0 1 0
             0 0 0 0 2 0 0 0 1 0 0
             0 0 0 0 1 0 0 1 0 0 0
             0 0 0 0 1 0 1 0 0 0 0
             0 0 0 0 0 1 0 0 1 0 0
             0 0 0 0 0 1 0 2 0 0 0
             0 0 0 0 0 1 2 0 0 0 0
             0 1 0 0 0 0 1 0 0 0 0
             0 2 0 0 0 0 0 0 1 0 0
             0 1 0 0 0 0 0 1 0 0 0
             2 0 0 0 0 0 0 0 1 0 0
             1 0 0 0 0 0 0 1 0 0 0
             1 0 0 0 0 0 1 0 0 0 0
             1 0 0 0 0 0 0 0 0 1 0
             0 1 0 0 0 0 0 0 0 1 0
             3 0 0 0 0 0 0 0 0 2 0
             0 0 2 0 0 0 0 0 1 0 0
             0 0 1 0 0 0 0 0 0 1 0
             0 0 1 0 0 0 0 1 0 0 0
             0 0 1 0 0 0 1 0 0 0 0]

    ions = [:NH4, :Na, :H, :Ca, :K, :Mg, :Cl, :NO3, :SO4, :HSO4, :OH]
    saltnames = [:CaNO32, :CaCl2, :CaSO4, :KHSO4, :K2SO4, :KNO3, :KCl, :MgSO4,
        :MgNO32, :MgCl2, :NaCl, :Na2SO4, :NaNO3, :NH42SO4, :NH4NO3,
        :NH4Cl, :NH4HSO4, :NaHSO4, :NH43HSO42, :H2SO4, :HHSO4, :HNO3, :HCl]

    # Check to make sure the the salts matrix corresponds to the salt chemical formulas
    derived_saltnames = [Symbol([r[i] > 0 ? Symbol(ions[i], r[i] > 1 ? r[i] : Symbol("")) :
                                 Symbol("")
                                 for i in 1:length(ions)]...) for r in eachrow(salts)]

    @test derived_saltnames == saltnames

    @test rank(salts) == 10

    f = qr(salts', ColumnNorm())

    pivotsalts = saltnames[f.p[1:rank(salts)]]

    # CaSO4 and MgSO4 have identical column norms, so QR pivoting may non-deterministically
    # choose either one. Both are valid pivot choices mathematically.
    expected_pivotsalts_base = [:NH43HSO42, :CaCl2, :K2SO4, :MgNO32, :Na2SO4,
        :H2SO4, :NH42SO4, :MgCl2, :NH4HSO4]
    @test length(pivotsalts) == 10
    @test issubset(expected_pivotsalts_base, pivotsalts)
    @test :MgSO4 in pivotsalts || :CaSO4 in pivotsalts

    nonpivotsalts = setdiff(saltnames, pivotsalts)
    expected_nonpivotsalts_base = [:CaNO32, :KHSO4, :KNO3, :KCl, :NaCl, :NaNO3,
        :NH4NO3, :NH4Cl, :NaHSO4, :HHSO4, :HNO3, :HCl]
    @test length(nonpivotsalts) == 13
    @test issubset(expected_nonpivotsalts_base, nonpivotsalts)
    # Exactly one of MgSO4 or CaSO4 should be in nonpivotsalts (the other is in pivotsalts)
    @test (:MgSO4 in nonpivotsalts) ⊻ (:CaSO4 in nonpivotsalts)

    @test rank(salts[f.p[1:rank(salts)], :]) == 10

    # Check how removing H2SO4, HHSO4, and CaSO4 affects the pivot salts,
    # as we don't use those ones in the reactions.
    saltnames_noh2so4 = saltnames[[1:2; 4:19; 22:23]]
    salts_noh2so4 = salts[[1:2; 4:19; 22:23], :]
    rank(salts_noh2so4)
    f = qr(salts_noh2so4', ColumnNorm())
    pivotsalts_noh2so4 = saltnames_noh2so4[f.p[1:rank(salts_noh2so4)]]
    # When CaSO4 is removed from the set, HNO3 becomes a pivot.
    # If MgSO4 was not a pivot before (CaSO4 was), then MgSO4 will also become a new pivot.
    diff_salts = setdiff(pivotsalts_noh2so4, pivotsalts)
    @test :HNO3 in diff_salts
    @test issubset(diff_salts, [:HNO3, :MgSO4])
end

@mtkmodel AqueousTestMolality begin
    @components begin
        aq = Aqueous()
    end
    @constants begin
        no_change = 0.0, [unit = u"mol/kg/s"]
        M_no_change = 0.0, [unit = u"mol/m^3/s"]
        m_one = 1.0, [unit = u"mol/kg"]
    end
    @equations begin
        # Fix the molalities at the initial concentration.
        D(aq.NH4.m) ~ 0
        D(aq.Na.m) ~ 0
        D(aq.H.m) ~ 0
        D(aq.Ca.m) ~ 0
        D(aq.K.m) ~ 0
        D(aq.Mg.m) ~ 0
        D(aq.Cl.m) ~ 0
        D(aq.NO3.m) ~ 0
        D(aq.SO4.m) ~ 0
        D(aq.HSO4.m) ~ 0

        # These species are neutral and are not in any salts, so
        # we fix their concentration.
        aq.NH3.m ~ 0
        aq.HCl_aq.m ~ 0
        aq.HNO3_aq.m ~ 0

        # We need to fix the molar concentration of one of the salts in
        # order to be able to calculate the water concentration.
        D(aq.NH42SO4.M) ~ M_no_change
    end
end

@named aqt = AqueousTestMolality()
sys = mtkcompile(aqt)

prob = ODEProblem(sys,
    [
        sys.aq.NH4.m => 1.0,
        sys.aq.Na.m => 1.0,
        sys.aq.H.m => 1.0,
        sys.aq.Ca.m => 0.5,
        sys.aq.K.m => 1.0,
        sys.aq.Mg.m => 0.5,
        sys.aq.Cl.m => 1.0,
        sys.aq.NO3.m => 1.0,
        sys.aq.SO4.m => 0.5,
        sys.aq.HSO4.m => 1.0, sys.aq.NH42SO4.M => 1.0e-8
    ],
    (0.0, 1.0), use_scc = false
)
sol = solve(prob, Rosenbrock23())

let
    vars = [
        sys.aq.NH4NO3.loga, sys.aq.NH4Cl.loga, sys.aq.NH4HSO4.loga, sys.aq.NH42SO4.loga,
        sys.aq.NH43HSO42.loga,
        sys.aq.NaCl.loga, sys.aq.Na2SO4.loga, sys.aq.NaNO3.loga, sys.aq.NaHSO4.loga,
        sys.aq.H2SO4.loga, sys.aq.HCl.loga, sys.aq.HNO3.loga, sys.aq.KHSO4.loga, sys.aq.HHSO4.loga,
        sys.aq.CaNO32.loga, sys.aq.CaCl2.loga, sys.aq.CaSO4.loga,
        sys.aq.KHSO4.loga, sys.aq.K2SO4.loga, sys.aq.KNO3.loga, sys.aq.KCl.loga,
        sys.aq.MgSO4.loga, sys.aq.MgNO32.loga, sys.aq.MgCl2.loga]
    collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

lb, ub = -20, 100

@mtkmodel AqueousTestActivity begin
    @components begin
        aq = Aqueous()
    end
    @constants begin
        M_no_change = 0.0, [unit = u"mol/m^3/s"]
        no_change = 0.0, [unit = u"mol/m^3/s"]
        m_no_change = 0.0, [unit = u"mol/kg/s"]
        m_one = 1.0, [unit = u"mol/kg"]
        a_one = 1.0, [unit = u"mol^2/kg^2"]
        a_one3 = 1.0, [unit = u"mol^3/kg^3"]
        m2_no_change = 0.0, [unit = u"(mol/kg)^2/s"]
        m3_no_change = 0.0, [unit = u"(mol/kg)^3/s"]
        m5_no_change = 0.0, [unit = u"(mol/kg)^5/s"]
    end
    @equations begin
        # aqueous activities set by Table 2.
        aq.CaNO32.loga_eq ~ clamp(13.0, lb, ub) # r1
        aq.CaCl2.loga_eq ~ clamp(27.0, lb, ub) # r2
        aq.CaSO4.loga_eq ~ clamp(-10.0, lb, ub) # r3
        aq.K2SO4.loga_eq ~ clamp(-4.2, lb, ub) # r4
        aq.KHSO4.loga_eq ~ clamp(3.2, lb, ub) # r5
        aq.KNO3.loga_eq ~ clamp(-0.14, lb, ub) # r6
        aq.KCl.loga_eq ~ clamp(2.2, lb, ub) # r7
        aq.MgSO4.loga_eq ~ clamp(12.0, lb, ub) # r8
        aq.MgNO32.loga_eq ~ clamp(35.0, lb, ub) # r9
        aq.MgCl2.loga_eq ~ clamp(51.0, lb, ub) # r10
        #aq.HNO3.loga_eq ~ clamp(-11.0, lb, ub) # r14
        #aq.HCl.loga_eq ~ clamp(14.0, lb, ub) # r16
        aq.Na2SO4.loga_eq ~ clamp(-0.73, lb, ub) # r19
        aq.NH42SO4.loga_eq ~ clamp(0.63, lb, ub) # r20
        aq.NaNO3.loga_eq ~ clamp(2.5, lb, ub) # r22
        aq.NaCl.loga_eq ~ clamp(3.6, lb, ub) # r23
        aq.NaHSO4.loga_eq ~ clamp(10.0, lb, ub) # r24
        aq.NH4HSO4.loga_eq ~ clamp(0.32, lb, ub) # r26
        aq.NH43HSO42.loga_eq ~ clamp(3.4, lb, ub) # r27

        D(aq.NH3_dissociated.M_eq) ~ 0

        # These species are neutral and are not in any salts, so
        # we fix their concentration.
        aq.NH3.m_eq ~ m_one
        aq.HCl_aq.m_eq ~ m_one
        aq.HNO3_aq.m_eq ~ m_one
        aq.H2O_dissociated.m_eq ~ m_one
    end
end

@named aqt = AqueousTestActivity()
sys = mtkcompile(aqt)

prob = ODEProblem(sys,
    [
        sys.aq.NH3_dissociated.M_eq => 1.0e-8
    ],
    initializealg = BrownFullBasicInit(nlsolve = RobustMultiNewton()),
    (0.0, 1.0), use_scc = false
)

sol = solve(prob, Rosenbrock23())

collect(zip(unknowns(sys), prob.u0))
collect(zip(unknowns(sys), sol.u[1]))

collect(zip([sys.aq.I, sys.aq.W], sol[[sys.aq.I, sys.aq.W]][1]))

salts = [:CaNO32, :CaCl2, :CaSO4, :KHSO4, :K2SO4, :KNO3, :KCl, :MgSO4, :MgNO32, :MgCl2,
    :NaCl, :Na2SO4, :NaNO3,
    :NH42SO4, :NH4NO3, :NH4Cl, :NH4HSO4, :NaHSO4, :NH43HSO42, :HHSO4, :HNO3, :HCl]

let
    vars = [reduce(getproperty, [sys, :aq, salt, :loga]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :m_eq]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :m_aq]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [sys.aq.NH4.m, sys.aq.Na.m, sys.aq.H.m, sys.aq.Ca.m, sys.aq.K.m, sys.aq.Mg.m,
        sys.aq.Cl.m, sys.aq.NO3.m, sys.aq.SO4.m, sys.aq.HSO4.m, sys.aq.OH.m, sys.aq.NH3.m,
        sys.aq.HNO3_aq.m, sys.aq.HCl_aq.m]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [sys.aq.NH4NO3.M, sys.aq.NH4Cl.M, sys.aq.NH4HSO4.M,
        sys.aq.NH42SO4.M, sys.aq.NH43HSO42.M,
        sys.aq.NaCl.M, sys.aq.Na2SO4.M, sys.aq.NaNO3.M, sys.aq.NaHSO4.M,
        sys.aq.HCl.M, sys.aq.HNO3.M, sys.aq.KHSO4.M, sys.aq.HHSO4.M,
        sys.aq.CaNO32.M, sys.aq.CaCl2.M, sys.aq.CaSO4.M,
        sys.aq.KHSO4.M, sys.aq.K2SO4.M, sys.aq.KNO3.M, sys.aq.KCl.M,
        sys.aq.MgSO4.M, sys.aq.MgNO32.M, sys.aq.MgCl2.M]
    collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [sys.aq.γ_NH4NO3, sys.aq.γ_NH4Cl, sys.aq.γ_NH4HSO4,
        sys.aq.γ_NH42SO4, sys.aq.γ_NH43HSO42,
        sys.aq.γ_NaCl, sys.aq.γ_Na2SO4, sys.aq.γ_NaNO3, sys.aq.γ_NaHSO4,
        sys.aq.γ_H2SO4, sys.aq.γ_HCl, sys.aq.γ_HNO3, sys.aq.γ_KHSO4, sys.aq.γ_HHSO4,
        sys.aq.γ_CaNO32, sys.aq.γ_CaCl2, sys.aq.γ_CaSO4,
        sys.aq.γ_KHSO4, sys.aq.γ_K2SO4, sys.aq.γ_KNO3, sys.aq.γ_KCl,
        sys.aq.γ_MgSO4, sys.aq.γ_MgNO32, sys.aq.γ_MgCl2]
    collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [sys.aq.maw_CaCl2.m_aw,
        sys.aq.maw_KHSO4.m_aw, sys.aq.maw_K2SO4.m_aw, sys.aq.maw_KNO3.m_aw, sys.aq.maw_KCl.m_aw,
        sys.aq.maw_MgSO4.m_aw, sys.aq.maw_MgNO32.m_aw, sys.aq.maw_MgCl2.m_aw,
        sys.aq.maw_NaCl.m_aw, sys.aq.maw_Na2SO4.m_aw, sys.aq.maw_NaNO3.m_aw,
        sys.aq.maw_NH42SO4.m_aw, sys.aq.maw_NH4NO3.m_aw, sys.aq.maw_NH4Cl.m_aw,
        sys.aq.maw_NH4HSO4.m_aw, sys.aq.maw_NH43HSO42.m_aw]
    collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [sys.aq.a_NH4NO3, sys.aq.a_NH4Cl, sys.aq.a_NH4HSO4,
        sys.aq.a_NH42SO4, sys.aq.a_NH43HSO42,
        sys.aq.a_NaCl, sys.aq.a_Na2SO4, sys.aq.a_NaNO3, sys.aq.a_NaHSO4,
        sys.aq.a_H2SO4, sys.aq.a_HCl, sys.aq.a_HNO3, sys.aq.a_KHSO4, sys.aq.a_HHSO4,
        sys.aq.a_CaNO32, sys.aq.a_CaCl2, sys.aq.a_CaSO4,
        sys.aq.a_KHSO4, sys.aq.a_K2SO4, sys.aq.a_KNO3, sys.aq.a_KCl,
        sys.aq.a_MgSO4, sys.aq.a_MgNO32, sys.aq.a_MgCl2]
    collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

lb = 1.0e-2
ub = 1.2

prob = ODEProblem(sys,
    [
        sys.aq.NH43HSO42.loga => log(30.0),
        sys.aq.CaCl2.loga => log(8.0e11),
        sys.aq.K2SO4.loga => log(0.016),
        sys.aq.MgNO32.loga => log(2.5e15),
        sys.aq.Na2SO4.loga => log(0.48),
        sys.aq.HNO3.loga => log(2.5e6),
        sys.aq.NH42SO4.loga => log(1.9),
        sys.aq.MgCl2.loga => log(9.6e21),
        sys.aq.MgSO4.loga => log(110000.0),
        sys.aq.NH4HSO4.loga => log(1.4), sys.aq.NH42SO4.M => 1.0e-8
    ],
    (0.0, 1.0), use_scc = false,
    guesses = [
        sys.aq.Ca.m => 0.5,
        sys.aq.SO4.m => 0.5,
        sys.aq.Mg.m => 0.5
    ]
)

sol = solve(prob, Rosenbrock23())

collect(zip(unknowns(sys), prob.u0))
collect(zip(unknowns(sys), sol.u[1]))

collect(zip([sys.aq.I, sys.aq.W], sol[[sys.aq.I, sys.aq.W]][1]))

unknowns(prob.f.initializeprob.f.sys)
equations(prob.f.initializeprob.f.sys)

iprob = prob.f.initializeprob
iprob.u0

ff(u, p) = iprob.f(abs.(u), p)
iiprob = NonlinearProblem(ff, iprob.u0, iprob.p)
sol = solve(iiprob, RobustMultiNewton())

guesses = unknowns(iprob.f.sys) .=> abs.(sol.u)

sol = solve(iprob, abstol = 5e-3, reltol = 5e-3)
sol.stats

collect(zip(
    unknowns(iprob.f.sys), round.(sol.u; sigdigits = 2), round.(sol.resid; sigdigits = 2)))

using SymbolicIndexingInterface: setp, getsym, parameter_values
using SciMLBase: remake

prob = remake(prob, u0 = [sys.aq.Ca.m => 1.0, sys.aq.CaCl2.M => 1.0])
f = getsym(prob, [sys.aq.γ_NaCl, sys.aq.γ_CaCl2, sys.aq.γ_NaNO3, sys.aq.γ_CaNO32])
f(prob)

f = getsym(prob,
    [sys.aq.Ca.M, sys.aq.Ca.m, sys.aq.CaCl2.M, sys.aq.W, sys.aq.maw_CaCl2.m_aw, sys.aq.I])
f(prob)

f = getsym(prob,
    [sys.aq.CaNO32.M, sys.aq.CaCl2.M, sys.aq.CaSO4.M, sys.aq.KHSO4.M,
        sys.aq.K2SO4.M, sys.aq.KNO3.M, sys.aq.KCl.M, sys.aq.MgSO4.M, sys.aq.MgNO32.M,
        sys.aq.MgCl2.M, sys.aq.NaCl.M, sys.aq.Na2SO4.M, sys.aq.NaNO3.M, sys.aq.NH42SO4.M,
        sys.aq.NH4NO3.M, sys.aq.NH4Cl.M, sys.aq.NH4HSO4.M, sys.aq.NH43HSO42.M, sys.aq.H2SO4.M,
        sys.aq.HHSO4.M, sys.aq.HNO3.M, sys.aq.HCl.M])
f(prob)

prob = remake(prob, u0 = [sys.aq.Cl.m => 1.0])
f = getsym(prob, [sys.aq.γ_NaCl, sys.aq.γ_CaCl2, sys.aq.γ_NaNO3, sys.aq.γ_CaNO32])
f(prob)

solve(prob, Rosenbrock23())

using ModelingToolkit
using ModelingToolkit: t_nounits, D_nounits
using OrdinaryDiffEq

@variables begin
    cation1m(t_nounits), [guess = 1.0]
    anion1m(t_nounits), [guess = 1.0]
    salt1M(t_nounits), [guess = 1.0]
    anion2m(t_nounits), [guess = 1.0]
    salt2M(t_nounits), [guess = 1.0]
    # a1(t_nounits), [guess=1.0]
    # a2(t_nounits), [guess=1.0]
    # I(t_nounits), [guess=1.0]
    W(t_nounits), [guess = 1.0]
    totalcat(t_nounits), [guess = 2.0]
    totalan1(t_nounits), [guess = 2.0]
    totalan2(t_nounits), [guess = 2.0]
    solid1(t_nounits), [guess = 0.0]
    solid2(t_nounits), [guess = 0.0]
end
@parameters begin
    maw1 = 2.0
    maw2 = 3.0
end
eqs = [
       # I ~ 0.5 * (cation1m + anion1m + cation1m + 4 * anion2m)
       # # a1 ~ cation1m * anion1m * I
       # a2 ~ cation1m * anion2m * I
       W ~ salt1M / maw1 + salt2M / maw2
       cation1m * W ~ salt1M + salt2M # Combining this equation with the next two makes, gives an identical equation to the charge balance equation below. This causes the solver to fail.
       anion1m * W ~ salt1M
       anion2m * W ~ salt2M
       totalcat ~ salt1M + salt2M
       totalan1 ~ salt1M + solid1 #+ 2solid2
       totalan2 ~ salt2M + 2solid1 #+ solid2
       #totalan2 ~ abs(totalan2)
       #0 ~ cation1m - anion1m - anion2m
       D_nounits(totalcat) ~ 0.0
       D_nounits(totalan1) ~ 0.0
       D_nounits(totalan2) ~ 0.0]
@named aqt_sys = System(eqs, t_nounits)
aqt = structural_simplify(aqt_sys)

prob = ODEProblem(aqt,
    [
        aqt.totalcat => 2.0,
        aqt.totalan1 => 1.5,
        aqt.totalan2 => 2.0
    ],
    (0.0, 1.0))

sol = solve(prob, Rosenbrock23())

sol[[aqt.salt1M, aqt.salt2M, aqt.W, aqt.solid1, aqt.solid2]]

unknowns(aqt)

using Plots
using Aerosol.ISORROPIA: Salt, Ion

@mtkmodel BinarySolution begin
    @components begin
        salt = Salt(drh = 0, l_t = 0, q, ν_cation, ν_anion,
            z_cation, z_anion)
    end
    @constants begin
        I_rate = 1.0, [unit = u"mol/kg/s"]
        W = 1.0e-6, [unit = u"kg/m^3"]
    end
    @variables begin
        Aᵧ_term(t), [description = "Debye-Hückel term used in Equation 7 and 8"]
        F_cat(t), [description = "Activity contribution from the cation"]
        F_an(t), [description = "Activity contribution from the anion"]
    end
    @equations begin
        D(salt.I) ~ I_rate
        salt.logm ~ log(salt.I / salt.m_one)
        salt.W ~ W

        Aᵧ_term ~ salt.Aᵧ * √salt.I / (√salt.I_one + √salt.I) / √salt.I_one
        F_cat ~ salt.Y * salt.logγ⁰ + Aᵧ_term * salt.zz * salt.Y
        F_an ~ salt.X * salt.logγ⁰ + Aᵧ_term * salt.zz * salt.X
        salt.F_cat ~ F_cat
        salt.F_an ~ F_an
    end
end

# Binary activity coefficients for comparison with Figures 1-4 in
# Kim, Y.P., Seinfeld, J.H. and Saxena, P., 1993. Atmospheric gas-aerosol equilibrium I.
# Thermodynamic model. Aerosol Science and Technology, 19(2), pp.157-181.

let
    salts = [
        (n = :NaCl, ν_cat = 1, ν_an = 1, z_cat = 1, z_an = abs(-1), q = 2.23),
        (n = :Na2SO4, ν_cat = 2, ν_an = 1, z_cat = 1, z_an = abs(-2), q = -0.19),
        (n = :NaNO3, ν_cat = 1, ν_an = 1, z_cat = 1, z_an = abs(-1), q = -0.39),
        (n = :NH4NO3, ν_cat = 1, ν_an = 1, z_cat = 1, z_an = abs(-1), q = -1.15)
    ]
    plots = []

    for s in salts
        @mtkcompile slt = BinarySolution(salt.z_cation = s.z_cat, salt.ν_cation = s.ν_cat,
            salt.z_anion = s.z_an, salt.ν_anion = s.ν_an, salt.q = s.q)

        prob = ODEProblem(slt,
            [
                slt.salt.I => 0.00001
            ],
            (0.0, 40.0), saveat = 0.1)
        sol = solve(prob, Rosenbrock23())

        p = plot(sol[slt.salt.I], exp.(sol[slt.salt.logγₜ₀]), label = :none,
            title = string(s.n), xlabel = "I (mol/kg)", ylabel = "γₜ₀",
            yscale = :log10)
        push!(plots, p)
        @info s.n, exp(sol[slt.salt.logγₜ₀][end])
    end
    plot(plots..., layout = (2, 2))
end
