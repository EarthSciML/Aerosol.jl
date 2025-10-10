using Aerosol
import Aerosol.ISORROPIA: Aqueous
using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
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

    @test pivotsalts == [:NH43HSO42, :CaCl2, :K2SO4, :MgNO32, :Na2SO4,
        :H2SO4, :NH42SO4, :MgCl2, :MgSO4, :NH4HSO4]

    nonpivotsalts = setdiff(saltnames, pivotsalts)
    @test nonpivotsalts == [:CaNO32, :CaSO4, :KHSO4, :KNO3, :KCl, :NaCl, :NaNO3,
        :NH4NO3, :NH4Cl, :NaHSO4, :HHSO4, :HNO3, :HCl]

    @test rank(salts[f.p[1:rank(salts)], :]) == 10


    # Check how removing H2SO4, HHSO4, and CaSO4 affects the pivot salts,
    # as we don't use those ones in the reactions.
    saltnames_noh2so4 = saltnames[[1:2; 4:19; 22:23]]
    salts_noh2so4 = salts[[1:2; 4:19; 22:23], :]
    rank(salts_noh2so4)
        f = qr(salts_noh2so4', ColumnNorm())
    pivotsalts_noh2so4 = saltnames_noh2so4[f.p[1:rank(salts_noh2so4)]]
    @test setdiff(pivotsalts_noh2so4, pivotsalts) == [:HNO3]
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
        D(aq.NH4.m) ~ no_change
        D(aq.Na.m) ~ no_change
        D(aq.H.m) ~ no_change
        D(aq.Ca.m) ~ no_change
        D(aq.K.m) ~ no_change
        D(aq.Mg.m) ~ no_change
        D(aq.Cl.m) ~ no_change
        D(aq.NO3.m) ~ no_change
        D(aq.SO4.m) ~ no_change
        D(aq.HSO4.m) ~ no_change

        # These species are neutral and are not in any salts, so
        # we fix their concentration.
        aq.NH3.m ~ m_one
        aq.HCl_aq.m ~ m_one
        aq.HNO3_aq.m ~ m_one

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
        sys.aq.HSO4.m => 1.0,

        sys.aq.NH42SO4.M => 1.0e-8,
    ],
    (0.0, 1.0), use_scc=false,
    )
sol = solve(prob, Rosenbrock23())

let
vars = [sys.aq.a_NH4NO3, sys.aq.a_NH4Cl, sys.aq.a_NH4HSO4, sys.aq.a_NH42SO4, sys.aq.a_NH43HSO42,
    sys.aq.a_NaCl, sys.aq.a_Na2SO4, sys.aq.a_NaNO3, sys.aq.a_NaHSO4,
    sys.aq.a_H2SO4, sys.aq.a_HCl, sys.aq.a_HNO3, sys.aq.a_KHSO4, sys.aq.a_HHSO4,
    sys.aq.a_CaNO32, sys.aq.a_CaCl2, sys.aq.a_CaSO4,
    sys.aq.a_KHSO4, sys.aq.a_K2SO4, sys.aq.a_KNO3, sys.aq.a_KCl,
    sys.aq.a_MgSO4, sys.aq.a_MgNO32, sys.aq.a_MgCl2]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end



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
        # This time, we fix the activities at their initial values.
        D(aq.a_NH43HSO42) ~ m5_no_change
        D(aq.a_CaCl2) ~ m3_no_change
        D(aq.a_K2SO4) ~ m3_no_change
        D(aq.a_MgNO32) ~ m3_no_change
        D(aq.a_Na2SO4) ~ m3_no_change
        D(aq.a_HNO3) ~ m2_no_change
        D(aq.a_NH42SO4) ~ m3_no_change
        D(aq.a_MgCl2) ~ m3_no_change
        D(aq.a_MgSO4) ~ m2_no_change
        D(aq.a_NH4HSO4) ~ m2_no_change

        # These species are neutral and are not in any salts, so
        # we fix their concentration.
        aq.NH3.m ~ m_one
        aq.HCl_aq.m ~ m_one
        aq.HNO3_aq.m ~ m_one

        # We need to fix the molar concentration of one of the salts in
        # order to be able to calculate the water concentration.
        D(aq.NH42SO4.M) ~ M_no_change
    end
end

@named aqt = AqueousTestActivity()
sys = mtkcompile(aqt)

prob = ODEProblem(sys,
    [
        sys.aq.a_NH43HSO42 => 0.079,
        sys.aq.a_CaCl2 => 1.4,
        sys.aq.a_K2SO4 => 0.039,
        sys.aq.a_MgNO32 => 0.4,
        sys.aq.a_Na2SO4 => 0.053,
        sys.aq.a_HNO3 => 1.0,
        sys.aq.a_NH42SO4 => 0.043,
        sys.aq.a_MgCl2 => 2.4,
        sys.aq.a_MgSO4 => 0.042,
        sys.aq.a_NH4HSO4 => 0.83,

        sys.aq.NH42SO4.M => 1.0e-8,
    ],
    (0.0, 1.0), use_scc=false,
    # guesses = [
    #     sys.aq.Ca.m => 0.5,
    #     sys.aq.SO4.m => 0.5,
    #     sys.aq.Mg.m => 0.5,
    # ]
    )

sol = solve(prob, Rosenbrock23())

collect(zip(unknowns(sys), prob.u0))
collect(zip(unknowns(sys), sol.u[1]))

collect(zip([sys.aq.I, sys.aq.W], sol[[sys.aq.I, sys.aq.W]][1]))

let
vars = [sys.aq.NH4.m, sys.aq.Na.m, sys.aq.H.m, sys.aq.Ca.m, sys.aq.K.m, sys.aq.Mg.m,
    sys.aq.Cl.m, sys.aq.NO3.m, sys.aq.SO4.m, sys.aq.HSO4.m, sys.aq.OH.m, sys.aq.NH3.m,
    sys.aq.HNO3_aq.m, sys.aq.HCl_aq.m]
collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
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

lb = 1.0e-2
ub = 1.5

prob = ODEProblem(sys,
    [
        sys.aq.a_NH43HSO42 =>  clamp(30.0, lb, ub),
        sys.aq.a_CaCl2 => clamp(8.0e11, lb, ub),
        sys.aq.a_K2SO4 => clamp(0.016, lb, ub),
        sys.aq.a_MgNO32 => clamp(2.5e15, lb, ub),
        sys.aq.a_Na2SO4 => clamp(0.48, lb, ub),
        sys.aq.a_HNO3 => clamp(2.5e6, lb, ub),
        sys.aq.a_NH42SO4 => clamp(1.9, lb, ub),
        sys.aq.a_MgCl2 => clamp(9.6e21, lb, ub),
        sys.aq.a_MgSO4 => clamp(110000.0, lb, ub),
        sys.aq.a_NH4HSO4 => clamp(1.4, lb, ub),

        sys.aq.NH42SO4.M => 1.0e-8,
    ],
        (0.0, 1.0), use_scc=false,
            guesses = [
        sys.aq.Ca.m => 0.5,
        sys.aq.SO4.m => 0.5,
        sys.aq.Mg.m => 0.5,
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
sol = solve(iiprob)

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

using Plots
using Aerosol.ISORROPIA: Salt, Ion

@mtkmodel NaClActivity begin
    @components begin
        Na = Ion(z=1)
        Cl = Ion(z=abs(-1))
        NaCl = Salt(cation = Na, anion = Cl, drh = 0.7528, l_t = 25.0, q = 2.23)
    end
    @constants begin
        I_rate = 1.0, [unit = u"mol/kg/s"]
    end
    @equations begin
        D(NaCl.I) ~ I_rate
        Na.m ~ NaCl.I
        Cl.m ~ NaCl.I
    end
end

@mtkcompile nacl = NaClActivity()

prob = ODEProblem(nacl,
    [
        nacl.NaCl.I => 0.1,
    ],
    (0.0, 40.0), saveat=1.0)
sol = solve(prob, Rosenbrock23())

plot(sol[nacl.NaCl.I], exp.(sol[nacl.NaCl.logγ⁰]))


@mtkmodel BinarySalt begin
    @components begin
        cation = Ion(z)
        anion = Ion(z)
        salt = Salt(cation = cation, anion = anion, drh = 0, l_t = 0, q)
    end
    @structural_parameters begin
        ν_cat = 1
        ν_an = 1
    end
    @constants begin
        I_rate = 1.0, [unit = u"mol/kg/s"]
    end
    @variables begin
        Aᵧ_term(t), [description = "Debye-Hückel term used in Equation 7 and 8"]
        F_cat(t), [description = "Activity contribution from the cation"]
        F_an(t), [description = "Activity contribution from the anion"]
    end
    @equations begin
        D(salt.I) ~ I_rate
        cation.m ~ ν_cat * salt.I
        anion.m ~ ν_an * salt.I

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

salts = [
    (n = :NaCl, ν_cat= 1, ν_an = 1, z_cat =1, z_an = abs(-1), q = 2.23),
    (n = :Na2SO4, ν_cat= 2, ν_an = 1, z_cat = 1, z_an = abs(-2), q = -0.19),
    (n = :NaNO3, ν_cat= 1, ν_an = 1, z_cat = 1, z_an = abs(-1), q = -0.39),
    (n = :NH4NO3, ν_cat= 1, ν_an = 1, z_cat = 1, z_an = abs(-1), q =  -1.15),
]
plots = []

for s in salts
    @mtkcompile slt = BinarySalt(cation.z=s.z_cat, ν_cat=s.ν_cat,
    anion.z=s.z_an, ν_an=s.ν_an, salt.q=s.q)

prob = ODEProblem(slt,
    [
        slt.salt.I => 0.00001,
    ],
    (0.0, 40.0), saveat=1.0)
sol = solve(prob, Rosenbrock23())

p = plot(sol[slt.salt.I], exp.(sol[slt.salt.logγₜ₀]), label=:none,
    title = string(s.n), xlabel = "I (mol/kg)", ylabel = "γₜ₀",
    yscale=:log10)
push!(plots, p)
end
plot(plots..., layout = (2,2))
