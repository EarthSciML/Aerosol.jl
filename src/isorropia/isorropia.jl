#module ISORROPIA
using EarthSciMLBase
using ModelingToolkit, Catalyst, Unitful
using IfElse

using Latexify
using NonlinearSolve
using Plots

export Isorropia

include("units.jl")
include("species.jl")
include("aqueous.jl")
include("deliquescence.jl")
include("solid.jl")
include("gas.jl")
include("water.jl")
include("reactions.jl")

"Function to transform variables to avoid negative concentrations."
trans(x) = abs(x)

"""
Substitute in transformed concentrations to avoid negative values, where 
`sys` is an equation system and `active_vars` is a vector of variables
to be substituted.
"""
function sub_trans(sys::ModelingToolkit.NonlinearSystem, active_vars::AbstractVector)
    trans_subs = [v => trans(v) for v ∈ active_vars]
    subbed_eqs = substitute(equations(sys), Dict(trans_subs))
    NonlinearSystem(subbed_eqs, states(sys), parameters(sys); name=nameof(sys))
end

"""
    Isorropia(rxn_nums)

An implementation of ISORROPIA II, a model for the thermodynamic equilibrium of gas-aerosol interactions, as described in:

> Fountoukis, C. and Nenes, A., 2007. ISORROPIA II: a computationally efficient thermodynamic equilibrium model for K+–Ca 2+–Mg 2+–NH 4+–Na+–SO 4 2−–NO 3−–Cl−–H 2 O aerosols. Atmospheric Chemistry and Physics, 7(17), pp.4639-4659.

"""
struct Isorropia <: EarthSciMLODESystem
    sys::NonlinearSystem
    Isorropia(sys::ModelingToolkit.NonlinearSystem) = new(sys)
    function Isorropia(rxn_nums::AbstractVector, name=:isorropia)
        rxns = reactions[rxn_nums]
        active_specs = active_species(rxns)
        active_vars = unique(vcat(vars.(active_specs)...))
        active_ions = intersect(active_specs, all_Ions)
        active_salts = intersect(active_specs, all_salts)
        active_solids = intersect(active_specs, all_Solids)
        active_gases = intersect(active_specs, all_Gases)

        water = Water(active_salts)
        ionic_strength = IonicStrength(active_ions, water.W)
        salt_act = Activities(active_salts, (s) -> activity(s, active_salts, ionic_strength.I, water.W))
        ion_act = Activities(active_ions, (i) -> activity(i, water.W))
        gas_act = Activities(active_gases, activity)
        solid_act = Activities(active_solids, activity)
        water_act = Activities([H2Oaq], activity)
        activities = extend(water_act, extend(salt_act, extend(ion_act, extend(gas_act, solid_act))))
        balance = Balance(T, active_vars, active_ions, active_salts)
        #out = Output(active_vars)

        sys = NonlinearSystem([], active_vars, [T, RH],
            systems=[rxn_sys(x, activities) for x ∈ rxns]; name=name,
        )
        sys = extend(sys, extend(ionic_strength, extend(activities, extend(balance, water))))
        sys = sub_trans(sys, active_vars)
        new(sys)
    end
end

"""
Create a system of equations for mass and charge balances.
"""
function Balance(T, active_vars, active_ions, active_salts)
    # Return the given species, or zero if the species is not in the list.
    g(v) = !isnothing(findfirst(isequal(v), active_vars)) ? v : 0
    all_active_ions = unique(vcat(active_ions, [[s.cation, s.anion] for s ∈ active_salts]...))

    @constants R = 8.31446261815324 [unit = u"m_air^3*Pa/K/mol", description = "Universal gas constant"]
    @constants PaPerAtm = 101325 [unit = u"Pa/atm", description = "Number of pascals per atmosphere"]
    atm2mol(p) = p * PaPerAtm / R / T

    totals = @parameters begin
        totalK = 1.0e-7, [unit = u"mol/m_air^3", description = "Total concentration of K"]
        totalCa = 1.0e-7, [unit = u"mol/m_air^3", description = "Total concentration of Ca"]
        totalMg = 1.0e-7, [unit = u"mol/m_air^3", description = "Total concentration of Mg"]
        totalNH = 1.0e-7, [unit = u"mol/m_air^3", description = "Total concentration of NH"]
        totalNa = 1.0e-7, [unit = u"mol/m_air^3", description = "Total concentration of Na"]
        totalSO4 = 1.0e-7, [unit = u"mol/m_air^3", description = "Total concentration of SO4"]
        totalNO3 = 1.0e-7, [unit = u"mol/m_air^3", description = "Total concentration of NO3"]
        totalCl = 1.0e-7, [unit = u"mol/m_air^3", description = "Total concentration of Cl"]
    end

    ratios = @variables begin
        R_1 = 1.0, [description = "Total sulfate ratio from Section 3.1 of Fountoukis and Nenes (2007)"]
        R_2 = 1.0, [description = "Crustal species and sodium ratio from Section 3.1 of Fountoukis and Nenes (2007)"]
        R_3 = 1.0, [description = "Crustal species ratio from Section 3.1 of Fountoukis and Nenes (2007)"]
    end

    eqs = [
        # Mass balances
        totalK ~ g(K_aq) + g(KHSO4_s) + 2g(K2SO4_s) + g(KNO3_s) + g(KCl_s)
        totalCa ~ g(Ca_aq) + g(CaSO4_s) + g(CaNO32_s) + g(CaCl2_s)
        totalMg ~ g(Mg_aq) + g(MgSO4_s) + g(MgNO32_s) + g(MgCl2_s)
        totalNH ~ g(NH4_aq) + g(NH3_aq) + atm2mol(g(NH3_g)) + g(NH4HSO4_s) +
                  2g(NH42SO4_s) + 3g(NH43HSO42_s) + g(NH4Cl_s) + g(NH4NO3_s)
        totalNa ~ g(Na_aq) + g(NaHSO4_s) + 2g(Na2SO4_s) + g(NaCl_s) + g(NaNO3_s)
        totalSO4 ~ g(SO4_aq) + g(HSO4_aq) +
                   g(KHSO4_s) + g(NaHSO4_s) + g(NH4HSO4_s) + g(CaSO4_s) + g(Na2SO4_s) + g(NH42SO4_s) +
                   2g(NH43HSO42_s) + g(K2SO4_s) + g(MgSO4_s)
        totalNO3 ~ g(NO3_aq) + g(HNO3_aq) + atm2mol(g(HNO3_g)) + g(NH4NO3_s) + g(NaNO3_s)
        totalCl ~ g(Cl_aq) + g(HCl_aq) + atm2mol(g(HCl_g)) +
                  g(NH4Cl_s) + g(NaCl_s) + 2g(CaCl2_s) + g(KCl_s) + 2g(MgCl2_s)
        # Charge balance
        0 ~ sum([i.m * i.z for i in all_active_ions])

        # Ratios from Section 3.1
        R_1 ~ (totalNH + totalCa + totalK + totalMg + totalNa) / totalSO4
        R_2 ~ (totalCa + totalK + totalMg + totalNa) / totalSO4
        R_3 ~ (totalCa + totalK + totalMg) / totalSO4
    ]
    nonzeros = [e.rhs !== 0 for e ∈ eqs]
    NonlinearSystem(eqs[nonzeros], ratios, totals, name=:balance)
end

"""
Create a system of equations for output variables.
"""
function Output(active_vars)
    names = Symbolics.tosymbol.(active_vars, escape=false)
    ovars = [only(@variables $(n) [unit=ModelingToolkit.get_unit(v), description="Output for $n"]) for 
            (n, v) ∈ zip(names, active_vars)]
    # This ends up with a double transform (e.g. ||v||), but is necessary to avoid cancelling.
    eqs = [ov ~ trans(v) for (ov, v) ∈ zip(ovars, active_vars)]
    NonlinearSystem(eqs, ovars, [], name=:out)
end

# Type 1 reactions (R_1 < 1.0)
rxn_nums = [3, 5, 11, 12, 13, 14, 15, 16, 17, 18, 24, 26]
sys = get_mtk(Isorropia(rxn_nums))

render(latexify(sys))
simple_sys = structural_simplify(sys)

render(latexify(simple_sys))
render(latexify(observed(simple_sys)))

ps = [simple_sys.totalSO4 => 1e-7, simple_sys.totalNa => 1e-20, simple_sys.totalNH => 1e-8,
    simple_sys.totalNO3 => 1e-7, simple_sys.totalCl => 1e-20, simple_sys.totalCa => 5e-9,
    simple_sys.totalK => 5e-9, simple_sys.totalMg => 1e-20]
prob = NonlinearProblem(simple_sys, [], ps)
sol = solve(prob, TrustRegion())
sol[sys.R_1] # Confirm that R_1 < 1.0
show(sol.resid) # Residuals should be small, but they are not. This means that the solver hasn't converged
sol[sys.W]
sol[sys.I]

sts = [sys.CaSO4_s, sys.Ca_aq, sys.SO4_aq, sys.KHSO4_s, sys.K_aq,
    sys.HSO4_aq, sys.H_aq, sys.NH3_g, sys.NH3_aq, sys.NH4_aq,
    sys.OH_aq, sys.HNO3_g, sys.NO3_aq, sys.HNO3_aq, sys.HCl_g,
    sys.Cl_aq, sys.HCl_aq, sys.NaHSO4_s, sys.Na_aq, sys.NH4HSO4_s]
show([v => trans(sol[v]) for v in sts]) # Negative numbers are fine because we take the absolute value.

sts = [sys.a_CaSO4_s, sys.a_KHSO4_s, sys.a_NaHSO4_s, sys.a_NH4HSO4_s, 
    sys.a_NH3_g, sys.a_HNO3_g, sys.a_HCl_g, sys.a_HSO4_aq, sys.a_H_aq, sys.a_SO4_aq, 
    sys.a_NH3_aq, sys.a_NH4_aq, sys.a_OH_aq, sys.a_HNO3_aq, sys.a_Cl_aq, sys.a_HCl_aq, 
    sys.a_CaSO4_aqs, sys.a_KHSO4_aqs, sys.a_HNO3_aqs, sys.a_NaHSO4_aqs, sys.a_NH4HSO4_aqs, sys.a_H2O]
show([v => sol[v] for v in sts])