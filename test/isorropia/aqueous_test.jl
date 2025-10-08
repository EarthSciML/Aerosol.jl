using Aerosol
import Aerosol.ISORROPIA: Aqueous
using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEq

#@named aq = Aqueous()

@mtkmodel AqueousTest begin
    @components begin
        aq = Aqueous()
    end
    @constants begin
        no_change = 0.0, [unit = u"mol/m^3/s"]
        m_no_change = 0.0, [unit = u"mol/kg/s"]
        m_one = 1.0, [unit = u"mol/kg"]
        a_one = 1.0, [unit = u"mol^2/kg^2"]
        a_one3 = 1.0, [unit = u"mol^3/kg^3"]
        m2_no_change = 0.0, [unit = u"(mol/kg)^2/s"]
        m3_no_change = 0.0, [unit = u"(mol/kg)^3/s"]
    end
    @variables begin
        #! format: off
        # TotalNH4(t), [unit = u"mol/m^3", description = "Total NH4 Molarity", guess = 1e-8]
        # TotalNa(t), [unit = u"mol/m^3", description = "Total Na+ Molarity", guess = 1e-8]
        # TotalCa(t), [unit = u"mol/m^3", description = "Total Ca2+ Molarity", guess = 1e-8]
        # TotalK(t), [unit = u"mol/m^3", description = "Total K+ Molarity", guess = 1e-8]
        # TotalMg(t), [unit = u"mol/m^3", description = "Total Mg2+ Molarity", guess = 1e-8]
        # TotalCl(t), [unit = u"mol/m^3", description = "Total Cl- Molarity", guess = 1e-8]
        # TotalNO3(t), [unit = u"mol/m^3", description = "Total NO3- Molarity", guess = 1e-8]
        # TotalSO4(t), [unit = u"mol/m^3", description = "Total SO4 2- Molarity", guess = 1e-8]
        # TotalHSO4(t), [unit = u"mol/m^3", description = "Total HSO4- Molarity", guess = 1e-8]
        # #! format: on

        # extraNH4(t), [unit = u"mol/m^3", description = "Extra NH4 mass", guess = 0.0]
        # extraNa(t), [unit = u"mol/m^3", description = "Extra Na mass", guess = 0.0]
        # extraCa(t), [unit = u"mol/m^3", description = "Extra Ca mass", guess = 0.0]
        # extraK(t), [unit = u"mol/m^3", description = "Extra K mass", guess = 0.0]
        # extraMg(t), [unit = u"mol/m^3", description = "Extra Mg mass", guess = 0.0]
        # extraCl(t), [unit = u"mol/m^3", description = "Extra Cl mass", guess = 0.0]
        # extraNO3(t), [unit = u"mol/m^3", description = "Extra NO3 mass", guess = 0.0]
        # #extraSO4(t), [unit = u"mol/m^3", description = "Extra SO4 mass", guess = 0.0]
        # extraHSO4(t), [unit = u"mol/m^3", description = "Extra HSO4 mass", guess = 0.0]
    end
    @equations begin
        # TotalNH4 ~ sum([aq.NH4NO3.M, aq.NH4Cl.M, aq.NH4HSO4.M, 2aq.NH42SO4.M, 3aq.NH43HSO42.M]) #+ extraNH4
        # TotalNa ~ sum([aq.NaCl.M, 2aq.Na2SO4.M, aq.NaNO3.M, aq.NaHSO4.M]) + extraNa
        # TotalCa ~ sum([aq.CaNO32.M, aq.CaCl2.M, aq.CaSO4.M]) + extraCa
        # TotalK ~ sum([aq.KHSO4.M, 2aq.K2SO4.M, aq.KNO3.M, aq.KCl.M]) + extraK
        # TotalMg ~ sum([aq.MgSO4.M, aq.MgNO32.M, aq.MgCl2.M]) + extraMg
        # TotalCl ~ sum([aq.NaCl.M, aq.KCl.M, 2aq.MgCl2.M, 2aq.CaCl2.M, aq.NH4Cl.M, aq.HCl.M]) #+ extraCl
        # TotalNO3 ~ sum([aq.NaNO3.M, aq.KNO3.M, 2aq.MgNO32.M, 2aq.CaNO32.M, aq.NH4NO3.M, aq.HNO3.M]) + extraNO3
        # TotalSO4 ~ sum([aq.Na2SO4.M, aq.K2SO4.M, aq.MgSO4.M, aq.CaSO4.M, aq.NH42SO4.M, aq.H2SO4.M,
        #     aq.NH43HSO42.M]) #+ extraSO4
        # TotalHSO4 ~ sum([aq.KHSO4.M, aq.NaHSO4.M, aq.NH4HSO4.M, aq.NH43HSO42.M, aq.HHSO4.M]) + extraHSO4


        # D(TotalNH4) ~ no_change
        # D(TotalNa) ~ no_change
        # D(TotalCa) ~ no_change
        # D(TotalK) ~ no_change
        # D(TotalMg) ~ no_change
        # D(TotalCl) ~ no_change
        # D(TotalNO3) ~ no_change
        # D(TotalSO4) ~ no_change
        # D(TotalHSO4) ~ no_change

        # aq.Cl.m ~ abs(aq.Cl.m)
        # aq.SO4.m ~ abs(aq.SO4.m)
        # aq.Na.m ~ abs(aq.Na.m)
        # aq.NH4.m ~ abs(aq.NH4.m)
        # aq.K.m ~ abs(aq.K.m)
        # aq.Mg.m ~ abs(aq.Mg.m)
        # aq.NO3.m ~ abs(aq.NO3.m)
        # aq.HSO4.m ~ abs(aq.HSO4.m)
        # aq.NH3.m ~ abs(aq.NH3.m)
        # aq.HNO3_aq.m ~ abs(aq.HNO3_aq.m)
        # aq.HCl_aq.m ~ abs(aq.HCl_aq.m)

        # Relationships between salts, just for testing purposes.
        # These are replaced by equilibrium equations in the real model.
       # aq.a_CaNO32 ~ a_one3 * 1e-8
       # aq.a_CaCl2 ~ a_one3 * 1e-8
       # aq.a_KHSO4 ~ a_one * 1e-8^2
        #aq.a_KNO3 ~ a_one * 1e-8
      #  aq.a_K2SO4 ~ a_one3 * 1e-8

       #aq.CaSO4.M ~ abs(aq.CaSO4.M)
        # aq.HCl.M ~ abs(aq.HCl.M)
        # aq.MgNO32.M ~ abs(aq.MgNO32.M)
        # aq.NH42SO4.M ~ abs(aq.NH42SO4.M)
        # aq.CaNO32.M ~ abs(aq.CaNO32.M)
        # aq.NH4Cl.M ~ abs(aq.NH4Cl.M)
        # aq.OH.m ~ abs(aq.OH.m)
        # aq.KNO3.M ~ abs(aq.KNO3.M)
        # aq.KCl.M ~ abs(aq.KCl.M)
        # aq.HNO3.M ~ abs(aq.HNO3.M)
        # aq.NaHSO4.M ~ abs(aq.NaHSO4.M)
        # aq.MgSO4.M ~ abs(aq.MgSO4.M)
        # aq.HHSO4.M ~ abs(aq.HHSO4.M)
        # aq.CaCl2.M ~ abs(aq.CaCl2.M)
        # aq.K2SO4.M ~ abs(aq.K2SO4.M)
        # aq.Na2SO4.M ~ abs(aq.Na2SO4.M)
        # aq.NaCl.M ~ abs(aq.NaCl.M)
        # aq.NH4NO3.M ~ abs(aq.NH4NO3.M)
        # aq.NH4HSO4.M ~ abs(aq.NH4HSO4.M)
        # aq.KHSO4.M ~ abs(aq.KHSO4.M)

        # aq.H.m ~ abs(aq.H.m)
        # aq.Na.m ~ abs(aq.Na.m)
        # aq.K.m ~ abs(aq.K.m)
        # aq.SO4.m ~ abs(aq.SO4.m)
        # aq.Ca.m ~ abs(aq.Ca.m)
        # aq.HSO4.m ~ abs(aq.HSO4.m)
        # aq.Mg.m ~ abs(aq.Mg.m)
        # aq.Cl.m ~ abs(aq.Cl.m)
        # aq.NO3.m ~ abs(aq.NO3.m)

        D(aq.a_CaNO32) ~ m3_no_change
        D(aq.a_CaCl2) ~ m3_no_change
        D(aq.a_KHSO4) ~ m2_no_change
        D(aq.a_K2SO4) ~ m3_no_change
        D(aq.a_KNO3) ~ m2_no_change
        D(aq.a_MgSO4) ~ m2_no_change
        D(aq.a_NaCl) ~ m2_no_change
        D(aq.a_NH42SO4) ~ m3_no_change
        D(aq.a_HNO3) ~ m2_no_change
        D(aq.a_NaHSO4) ~ m2_no_change

        #aq.OH.m ~ m_one

        aq.NH3.m ~ m_one
        aq.HCl_aq.m ~ m_one
        aq.HNO3_aq.m ~ m_one
    end
end

@named aqt = AqueousTest()
sys = mtkcompile(aqt)
sys2 = mtkcompile(aqt, fully_determined = false)
equations(sys)
unknowns(sys)

prob = ODEProblem(sys,
    [
#        sys.aq.T => 293.15, sys.aq.RH => 0.3,

        # sys.aq.a_CaNO32 ~ 1.0e-16,
        # sys.aq.a_CaCl2 ~ 1.0e-16,
        # sys.aq.a_KHSO4 ~ 1.0e-16,
        # sys.aq.a_K2SO4 ~ 1.0e-16,
        # sys.aq.a_KNO3 ~ 1.0e-16,
        # sys.aq.a_MgSO4 ~ 1.0e-16,
        # sys.aq.a_NaCl ~ 1.0e-16,
        # sys.aq.a_NH42SO4 ~ 1.0e-16,
        # sys.aq.a_HNO3 ~ 1.0e-16,
        # sys.aq.a_NaHSO4 ~ 1.0e-16,

        #  sys.aq.NH4.m ~ 1.0e-8,
        # sys.aq.NO3.m ~ 1.0e-8,
        # sys.aq.H.m ~ 1.0e-8,
        # sys.aq.Na.m ~ 1.0e-8,
        # sys.aq.K.m ~ 1.0e-8,
        # sys.aq.SO4.m ~ 1.0e-8,
        # sys.aq.Ca.m ~ 1.0e-8,
        # sys.aq.HSO4.m ~ 1.0e-8,
        # sys.aq.Mg.m ~ 1.0e-8,
        # sys.aq.Cl.m ~ 1.0e-8,

        # sys.TotalNH4 => 1e-8,
        # sys.TotalNa => 1e-8,
        # sys.TotalCl => 1e-8,
        # sys.TotalSO4 => 1e-8,
        # sys.TotalK => 1e-8,
        # sys.TotalMg => 1e-8,
        # sys.TotalNO3 => 1e-8,
        # sys.TotalHSO4 => 1e-8,
        # sys.TotalCa => 1e-8,

        #sys.aq.H.m => 1.0,
    #    sys.aq.Ca.m => 1.0,
      #  sys.aq.NO3.m => 1.0,
      #  sys.aq.HSO4.m => 1.0,
      #  sys.aq.NH3.m => 1.0,
        # sys.aq.HNO3_aq.m => 1.0,
        # sys.aq.HCl_aq.m => 1.0
    ],
    (0.0, 1.0), use_scc = false,
    )
sol = solve(prob, abstol=5e-3, reltol = 5e-3)

collect(zip(unknowns(sys), prob.u0))
collect(zip(unknowns(sys), sol.u[1]))

collect(zip([sys.aq.I, sys.aq.W], sol[[sys.aq.I, sys.aq.W]][1]))

let
vars = [sys.aq.NH4.m, sys.aq.Na.m, sys.aq.H.m, sys.aq.Ca.m, sys.aq.K.m, sys.aq.Mg.m,
    sys.aq.Cl.m, sys.aq.NO3.m, sys.aq.SO4.m, sys.aq.HSO4.m, sys.aq.OH.m, sys.aq.NH3.m,
    sys.aq.HNO3_aq.m, sys.aq.HCl_aq.m]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
vars = [sys.aq.NH4NO3.M, sys.aq.NH4Cl.M, sys.aq.NH4HSO4.M, sys.aq.NH42SO4.M, sys.aq.NH43HSO42.M,
    sys.aq.NaCl.M, sys.aq.Na2SO4.M, sys.aq.NaNO3.M, sys.aq.NaHSO4.M,
    sys.aq.H2SO4.M, sys.aq.HCl.M, sys.aq.HNO3.M, sys.aq.KHSO4.M, sys.aq.HHSO4.M,
    sys.aq.CaNO32.M, sys.aq.CaCl2.M, sys.aq.CaSO4.M,
    sys.aq.KHSO4.M, sys.aq.K2SO4.M, sys.aq.KNO3.M, sys.aq.KCl.M,
    sys.aq.MgSO4.M, sys.aq.MgNO32.M, sys.aq.MgCl2.M]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
vars = [sys.aq.γ_NH4NO3, sys.aq.γ_NH4Cl, sys.aq.γ_NH4HSO4, sys.aq.γ_NH42SO4, sys.aq.γ_NH43HSO42,
    sys.aq.γ_NaCl, sys.aq.γ_Na2SO4, sys.aq.γ_NaNO3, sys.aq.γ_NaHSO4,
    sys.aq.γ_H2SO4, sys.aq.γ_HCl, sys.aq.γ_HNO3, sys.aq.γ_KHSO4, sys.aq.γ_HHSO4,
    sys.aq.γ_CaNO32, sys.aq.γ_CaCl2, sys.aq.γ_CaSO4,
    sys.aq.γ_KHSO4, sys.aq.γ_K2SO4, sys.aq.γ_KNO3, sys.aq.γ_KCl,
    sys.aq.γ_MgSO4, sys.aq.γ_MgNO32, sys.aq.γ_MgCl2]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
vars = [sys.aq.a_NH4NO3, sys.aq.a_NH4Cl, sys.aq.a_NH4HSO4, sys.aq.a_NH42SO4, sys.aq.a_NH43HSO42,
    sys.aq.a_NaCl, sys.aq.a_Na2SO4, sys.aq.a_NaNO3, sys.aq.a_NaHSO4,
    sys.aq.a_H2SO4, sys.aq.a_HCl, sys.aq.a_HNO3, sys.aq.a_KHSO4, sys.aq.a_HHSO4,
    sys.aq.a_CaNO32, sys.aq.a_CaCl2, sys.aq.a_CaSO4,
    sys.aq.a_KHSO4, sys.aq.a_K2SO4, sys.aq.a_KNO3, sys.aq.a_KCl,
    sys.aq.a_MgSO4, sys.aq.a_MgNO32, sys.aq.a_MgCl2]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end


unknowns(prob.f.initializeprob.f.sys)
equations(prob.f.initializeprob.f.sys)

iprob = prob.f.initializeprob
sol = solve(iprob, abstol=5e-3, reltol = 5e-3)
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
        cation1m(t_nounits), [guess=1.0]
        anion1m(t_nounits), [guess=1.0]
        salt1M(t_nounits), [guess=1.0]
        anion2m(t_nounits), [guess=1.0]
        salt2M(t_nounits), [guess=1.0]
        # a1(t_nounits), [guess=1.0]
        # a2(t_nounits), [guess=1.0]
        # I(t_nounits), [guess=1.0]
        W(t_nounits), [guess=1.0]
        totalcat(t_nounits), [guess=2.0]
        totalan1(t_nounits), [guess=2.0]
        totalan2(t_nounits), [guess=2.0]
        solid1(t_nounits), [guess=0.0]
        solid2(t_nounits), [guess=0.0]
    end
    @parameters begin
        maw1=2.0
        maw2=3.0
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
        D_nounits(totalan2) ~ 0.0
    ]
@named aqt_sys = System(eqs, t_nounits)
aqt = structural_simplify(aqt_sys)

prob = ODEProblem(aqt,
    [
        aqt.totalcat => 2.0,
        aqt.totalan1 => 1.5,
        aqt.totalan2 => 2.0,
    ],
    (0.0, 1.0))

sol = solve(prob, Rosenbrock23())

sol[[aqt.salt1M, aqt.salt2M, aqt.W, aqt.solid1, aqt.solid2]]

unknowns(aqt)
