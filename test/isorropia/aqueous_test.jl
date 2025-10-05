using Aerosol
import Aerosol.ISORROPIA: Aqueous
using ModelingToolkit
using ModelingToolkit: D
using DynamicQuantities
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve

@named aq = Aqueous()

@mtkmodel AqueousTest begin
    @components begin
        aq = Aqueous()
    end
    @constants begin
        no_change = 0.0, [unit = u"mol/m^3/s"]
        M_one = 1.0, [unit = u"mol/m^3"]
    end
    @equations begin
        D(aq.Ca.M) ~ no_change
        D(aq.Cl.M) ~ no_change
        D(aq.SO4.M) ~ no_change
        D(aq.Na.M) ~ no_change
        D(aq.NH4.M) ~ no_change
        D(aq.K.M) ~ no_change
        D(aq.Mg.M) ~ no_change
        D(aq.NO3.M) ~ no_change
        D(aq.HSO4.M) ~ no_change
        D(aq.NH3.M) ~ no_change
        D(aq.HNO3_aq.M) ~ no_change
        D(aq.HCl_aq.M) ~ no_change

        # Relationships between salts, just for testing purposes.
        # These are replaced by equilibrium equations in the real model.
        aq.CaSO4.M ~ aq.NH4HSO4.M
        aq.HCl.M ~ aq.NaCl.M
        aq.MgNO32.M ~ aq.Na2SO4.M
        aq.NH42SO4.M ~ aq.NH4NO3.M
        aq.CaNO32.M ~ aq.CaCl2.M
        aq.NH4Cl.M ~ aq.MgCl2.M
        aq.OH.M ~ 0
        aq.KNO3.M ~ aq.NH43HSO42.M
        aq.KCl.M + 2aq.HNO3.M ~ aq.K2SO4.M
        aq.HNO3.M ~ aq.NaNO3.M
        aq.NaHSO4.M ~ 0
        aq.MgSO4.M ~ 0
    end
end

@named aqt = AqueousTest()
sys = mtkcompile(aqt)
unknowns(sys)

prob = ODEProblem(sys, [sys.aq.T => 293.15, sys.aq.RH => 0.3], (0.0, 1.0))
prob = ODEProblem(sys, [
        sys.aq.T => 293.15, sys.aq.RH => 0.3,
        sys.aq.Ca.M => 1.0,
        sys.aq.Cl.M => 1.0,
        sys.aq.SO4.M => 1.0,
        sys.aq.Na.M => 1.0,
        sys.aq.NH4.M => 1.0,
        sys.aq.K.M => 1.0,
        sys.aq.Mg.M => 1.0,
        sys.aq.NO3.M => 1.0,
        sys.aq.HSO4.M => 1.0,
        sys.aq.NH3.M => 1.0,
        sys.aq.HNO3_aq.M => 1.0,
        sys.aq.HCl_aq.M => 1.0
    ], (0.0, 1.0), use_scc=false)

sol = solve(prob, Rosenbrock23())

sol[[sys.aq.I, sys.aq.W]]


unknowns(prob.f.initializeprob.f.sys)
equations(prob.f.initializeprob.f.sys)

isys = mtkcompile(ModelingToolkit.generate_initializesystem(sys))

equations(isys)

unknowns(isys)

defaults(sys)

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
