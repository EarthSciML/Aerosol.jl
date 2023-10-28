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

"""
    Isorropia(t)

An implementation of ISORROPIA II, a model for the thermodynamic equilibrium of gas-aerosol interactions, as described in:

> Fountoukis, C. and Nenes, A., 2007. ISORROPIA II: a computationally efficient thermodynamic equilibrium model for K+–Ca 2+–Mg 2+–NH 4+–Na+–SO 4 2−–NO 3−–Cl−–H 2 O aerosols. Atmospheric Chemistry and Physics, 7(17), pp.4639-4659.

"""
# struct Isorropia <: EarthSciMLODESystem
#     sys::ODESystem
#     rxn_sys::ReactionSystem
#     Isorropia(sys::ModelingToolkit.ODESystem, rxn_sys::ReactionSystem) = new(sys, rxn_sys)
#     function Isorropia(tt; name=:Isorropia)
#         global t = tt # Make t global so that it can be used when defining MTK variables.



# Load the other files into the global namespace. Ideally we would use functions 
# instead, but there is so much systemic complexity it's difficult to encapsulate
# everything into functions.

rxn_nums = [3, 5, 11, 12, 13, 14, 15, 16, 17, 18, 24, 26]
type1_rxns = reactions[rxn_nums]

active_specs = active_species(type1_rxns)
active_vars = unique(vcat(vars.(active_specs)...))
active_ions = intersect(active_specs, all_Ions)
active_salts = intersect(active_specs, all_salts)
active_solids = intersect(active_specs, all_Solids)
active_gases = intersect(active_specs, all_Gases)

ionic_strength = IonicStrength(active_ions)
water = Water(active_salts)
salt_act = SaltActivities(active_salts, ionic_strength.I)
ion_act = Activities(active_ions)
gas_act = Activities(active_gases)
solid_act = Activities(active_solids)
water_act = Activities([H2Oaq])
activities = extend(water_act, extend(salt_act, extend(ion_act, extend(gas_act, solid_act))))
balance = Balance(water.W, T, active_vars, active_ions)

length(type1_rxns) + length(equations(balance)) == length(active_vars)

@named sys = NonlinearSystem(
    [],
    active_vars,
    [T, RH],
    systems=[rxn_sys(x, activities) for x ∈ type1_rxns],
)
sys = extend(water, extend(activities, extend(balance, sys)))
render(latexify(sys))
sys = structural_simplify(sys)

render(latexify(sys))
render(latexify(observed(sys)))

prob = NonlinearProblem(sys, ModelingToolkit.get_defaults(sys), ModelingToolkit.get_defaults(sys))
sol = solve(prob, TrustRegion())
[u => sol[u] for u in states(sys)]
plot(sol)

@variables R_1 = 1.0 [description = "Total sulfate ratio from Section 3.1 of Fountoukis and Nenes (2007)"]
@variables R_2 = 1.0 [description = "Crustal species and sodium ratio from Section 3.1 of Fountoukis and Nenes (2007)"]
@variables R_3 = 1.0 [description = "Crustal species ratio from Section 3.1 of Fountoukis and Nenes (2007)"]

statevars = [all_solids; all_ions; all_gases; I; W]
ps = [T; RH; H2O_aq; metastable; totals]

#@parameters so4rate=100 [unit = u"s^-1", description = "Rate of SO4 to conversion aerosol phase (pseudo-instantaneous)"] # should be a @constants
eqs = vcat(
    # All SO4 immediately goes to aerosol phase as per Section 3.3 (item 1) of Fountoukis and Nenes (2007).
    0 ~ SO4_g / SO4_aq,

    # Ratios from Section 3.1
    R_1 ~ (totalNH + totalCa + totalK + totalMg + totalNa) / totalSO4,
    R_2 ~ (totalCa + totalK + totalMg + totalNa) / totalSO4,
    R_3 ~ (totalCa + totalK + totalMg) / totalSO4, all_rxn_eqs...)

@named ns = NonlinearSystem(eqs, statevars, ps, systems=[IW, DRH, [rxn.sys for rxn ∈ all_rxns]...])
structural_simplify(ns)

# rxn_sys = ReactionSystem(eqs, t, statevars, [so4rate; ps]; 
#     systems=[IW, DRH, [rxn.sys for rxn ∈ all_rxns]...], 
#     checks=true, combinatoric_ratelaws=false, name=name)
# sys = convert(ODESystem, rxn_sys)
# new(sys, rxn_sys)
#     end
# end
# end

"""
Create a system of equations for mass and charge balances.
"""
function Balance(W, T, active_vars, active_ions)
    # Return the given species, or zero if the species is not in the list.
    g(v) = !isnothing(findfirst(isequal(v), active_vars)) ? v : 0

    @constants R = 8.31446261815324 [unit = u"m_air^3*Pa/K/mol", description = "Universal gas constant"]
    @constants PaPerAtm = 101325 [unit = u"Pa/atm", description = "Number of pascals per atmosphere"]
    atm2mol(p) = p * PaPerAtm / R / T

    @parameters totalK = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of K"]
    @parameters totalCa = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of Ca"]
    @parameters totalMg = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of Mg"]
    @parameters totalNH = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of NH"]
    @parameters totalNa = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of Na"]
    @parameters totalSO4 = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of SO4"]
    @parameters totalNO3 = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of NO3"]
    @parameters totalCl = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of Cl"]
    totals = [totalK, totalCa, totalMg, totalNH, totalNa, totalSO4, totalNO3, totalCl]

    eqs = [
        # Mass balances
        totalK ~ g(K_aq) * W + g(KHSO4_s) + 2g(K2SO4_s) + g(KNO3_s) + g(KCl_s)
        totalCa ~ g(Ca_aq) * W + g(CaSO4_s) + g(CaNO32_s) + g(CaCl2_s)
        totalMg ~ g(Mg_aq) * W + g(MgSO4_s) + g(MgNO32_s) + g(MgCl2_s)
        totalNH ~ (g(NH4_aq) + g(NH3_aq)) * W + atm2mol(g(NH3_g)) + g(NH4HSO4_s) +
                  2g(NH42SO4_s) + 3g(NH43HSO42_s) + g(NH4Cl_s) + g(NH4NO3_s)
        totalNa ~ g(Na_aq) * W + g(NaHSO4_s) + 2g(Na2SO4_s) + g(NaCl_s) + g(NaNO3_s)
        totalSO4 ~ (g(SO4_aq) + g(HSO4_aq)) * W +
                g(KHSO4_s) + g(NaHSO4_s) + g(NH4HSO4_s) + g(CaSO4_s) + g(Na2SO4_s) + g(NH42SO4_s) +
                   2g(NH43HSO42_s) + g(K2SO4_s) + g(MgSO4_s)
        totalNO3 ~ (g(NO3_aq) + g(HNO3_aq)) * W + atm2mol(g(HNO3_g)) + g(NH4NO3_s) + g(NaNO3_s)
        totalCl ~ (g(Cl_aq) + g(HCl_aq)) * W + atm2mol(g(HCl_g)) +
                g(NH4Cl_s) + g(NaCl_s) + 2g(CaCl2_s) + g(KCl_s) + 2g(MgCl2_s)
        # Charge balance
        0 ~ sum([i.m * i.z for i in active_ions])
    ]
    nonzeros = [e.rhs !== 0 for e ∈ eqs]
    NonlinearSystem(eqs[nonzeros], [], totals, name=:balance)
end

CaSO4_s, Ca_aq, SO4_aq, KHSO4_s, K_aq, HSO4_aq, H_aq, NH3_g, 
NH3_aq, NH4_aq, OH_aq, HNO3_g, NO3_aq, 
HNO3_aq, HCl_g, Cl_aq, HCl_aq, NaHSO4_s, Na_aq, NH4HSO4_s, SO4_g




# TODO(CT): Missing H2O_g from table 3
kept_vars = [
    NaHSO4_s, NH4HSO4_s, KHSO4_s, CaSO4_s,
    Na_aq, NH4_aq, H_aq, HSO4_aq, SO4_aq, NO3_aq, Cl_aq, Ca_aq, K_aq, H2O_aq, # H2O_g
    NH3_g, NO3_aq, Cl_aq, NH3_aq, HNO3_aq, HCl_aq, # Minor species
]
push!(kept_vars, I, W)

sub_rules = [v => 0.0 for v ∈ setdiff(states(ns), kept_vars)]

eqS1 = [substitute(eq, Dict(sub_rules)) for eq ∈ eqs]

render(latexify(eqS1))

statevarsS1 = intersect(statevars, kept_vars)

function num_variables(eq)
    # Get the sources of each variable in the equation
    varsource = [getmetadata(var, Symbolics.VariableSource)[1] for var ∈ get_variables(eq)]
    # Count the number of state variables in the equation, where `varsource` is :variables and 
    # not :parameters or :constants.
    if length(varsource) == 0
        return 0
    end
    sum(varsource .== :variables)
end

eqs1_filtered = eqS1[[num_variables(eq) > 0 for eq ∈ eqS1]]

render(latexify(eqs1_filtered))

@named ns = NonlinearSystem(eqs1_filtered, statevarsS1, ps)
structural_simplify(ns)