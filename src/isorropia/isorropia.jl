#using Revise
using ModelingToolkit, NonlinearSolve, Unitful, Latexify
using IfElse
using Test

# Register the units we'll be using
module MyUnits
using Unitful
@unit m_air "m(air)" MAir 1u"m" false
@unit kg_water "kg(water)" KgWater 1u"kg" false
end
Unitful.register(MyUnits)

# Specify concentration bounds
ubound = 1.0e6
lbound = 1.0e-10

# Miscellaneous variables and parameters
@variables I = 1.0 [bounds=(lbound, ubound), unit = u"mol/kg_water", description = "Ionic strength"]
@variables W = 1.0e-5 [bounds=(lbound, ubound), unit = u"kg_water/m_air^3", description = "Aerosol water content"]

@parameters T = 293.15 [unit = u"K", description = "Temperature"]
@parameters RH = 0.3 [description = "Relative humidity (expressed on a scale from 0 to 1)"] # unitless

"""
A species represents a chemical species in the system.

Chemical species should have an `activity` function which returns the activity
of the species as specified in Section 2.2 of Fountoukis and Nenes (2007). 
"""
abstract type Species end
activity(s::Species) = error("activity not defined for $(typeof(s))")

""" 
The activity of multiple species is the product of their activities
as shown in Table 2 of Fountoukis and Nenes (2007).
"""
activity(s::AbstractVector) = reduce(*, activity.(s))

# Load the other files
include("isorropia_aqueous.jl")
include("isorropia_solid.jl")
include("isorropia_gas.jl")
include("isorropia_water.jl")

"""
Define a reaction based on information in Table 2 of Fountoukis and Nenes (2007).
"""
struct Reaction
    reactant
    product
    sys

    function Reaction(reactant, product, K⁰::Number, K⁰units, hgroup::Number, cgroup::Number; name="rxn")
        @constants K⁰ = K⁰ [unit = K⁰units, description = "Equilibrium constant at 298.15 K"]
        @constants H_group = hgroup [description = "ΔH⁰ / (R * T₀)"] # unitless
        @constants C_group = cgroup [description = "ΔC⁰ₚ / R"] # unitless
        sys = NonlinearSystem([], [], [T₀, K⁰, H_group, C_group]; name=name)
        new(reactant, product, sys)
    end
end

# Equation 5: Equilibrium constant
@constants T₀ = 293.15 [unit = u"K", description = "Standard temperature"]
K_eq(r::Reaction) = r.sys.K⁰ * exp(-r.sys.H_group * (r.sys.T₀ / T - 1) - r.sys.C_group * (1 + log(r.sys.T₀ / T) - r.sys.T₀ / T))

"""
Assemble an equation for a reaction based on Table 2 of Fountoukis and Nenes (2007), where
the left-hand side is the equilibrium constant and the right-hand side is the activity.
"""
function equation(r::Reaction)
    K_eq(r) ~ activity(r.product) / activity(r.reactant)
end

#== 
TODO: Assuming that H_group and C_group are zero when they are left out of Table 2. 
Is this correct?
==#
@named rxn1 = Reaction(CaNO32s, CaNO32_aqs, 6.067e5, u"mol^3/kg^3", -11.299, 0.0)
@named rxn2 = Reaction(CaCl2s, CaCl2_aqs, 7.974e11, u"mol^3/kg^3", -14.087, 0.0)
# Ignoring reaction 3 because activity coefficient of CaSO4 is zero.
#@named rxn3 = Reaction(CaSO4s, CaSO4_aqs, 4.319e-5, u"mol^2/kg^2", 0.0, 0.0)
@named rxn4 = Reaction(K2SO4s, K2SO4_aqs, 1.569e-2, u"mol^3/kg_water^3", -9.589, 45.807)
@named rxn5 = Reaction(KHSO4s, KHSO4_aqs, 24.016, u"mol^2/kg_water^2", -8.423, 17.964)
@named rxn6 = Reaction(KNO3s, KNO3_aqs, 0.872, u"mol^2/kg_water^2", 14.075, 19.388)
@named rxn7 = Reaction(KCls, KCl_aqs, 8.680, u"mol^2/kg_water^2", -6.167, 19.953)
@named rxn8 = Reaction(MgSO4s, MgSO4_aqs, 1.079e5, u"mol^2/kg_water^2", 36.798, 0.0)
@named rxn9 = Reaction(MgNO32s, MgNO32_aqs, 2.507e15, u"mol^3/kg_water^3", -8.754, 0.0)
@named rxn10 = Reaction(MgCl2s, MgCl2_aqs, 9.557e21, u"mol^3/kg_water^3", -1.347, 0.0)
@named rxn11 = Reaction(HSO4_ion, [H_ion, SO4_ion], 1.015e-2, u"mol/kg_water", 8.85, 25.14)
@named rxn12 = Reaction(NH3g, NH3_ion, 5.764e1, u"mol/kg_water/atm", 13.79, -5.39)
@named rxn13 = Reaction([NH3_ion, H2O_aq], [NH4_ion, OH_ion], 1.805e-5, u"mol/kg_water", -1.50, 26.92)
@named rxn14 = Reaction(HNO3g, HNO3_aqs, 2.511e6, u"mol^2/kg_water^2/atm", 29.17, 16.83)
@named rxn15 = Reaction(HNO3g, HNO3_ion, 2.1e5, u"mol/kg_water/atm", 29.17, 16.83)
@named rxn16 = Reaction(HClg, [H_ion, Cl_ion], 1.971e6, u"mol^2/kg_water^2/atm", 30.20, 19.91)
@named rxn17 = Reaction(HClg, HCl_ion, 2.5e3, u"mol/kg_water/atm", 30.20, 19.91)
@named rxn18 = Reaction(H2O_aq, [H_ion, OH_ion], 1.010e-14, u"mol^2/kg_water^2", -22.52, 26.92)
@named rxn19 = Reaction(Na2SO4s, Na2SO4_aqs, 4.799e-1, u"mol^3/kg_water^3", 0.98, 39.75)
@named rxn20 = Reaction(NH42SO4s, NH42SO4_aqs, 1.87e0, u"mol^3/kg_water^3", -2.65, 38.57)
@named rxn21 = Reaction(NH4Cls, [NH3g, HClg], 1.086e-16, u"atm^2", -71.00, 2.40)
@named rxn22 = Reaction(NaNO3s, NaNO3_aqs, 1.197e1, u"mol^2/kg_water^2", -8.22, 16.01)
@named rxn23 = Reaction(NaCls, NaCl_aqs, 3.766e1, u"mol^2/kg_water^2", -1.56, 16.90)
@named rxn24 = Reaction(NaHSO4s, NaHSO4_aqs, 2.413e4, u"mol^2/kg_water^2", 0.79, 14.75)
@named rxn25 = Reaction(NH4NO3s, [NH3g, HNO3g], 4.199e-17, u"atm^2", -74.375, 6.025)
@named rxn26 = Reaction(NH4HSO4s, NH4HSO4_aqs, 1.383e0, u"mol^2/kg_water^2", -2.87, 15.83)
@named rxn27 = Reaction(NH43HSO42s, NH43HSO42_aqs, 2.972e1, u"mol^5/kg_water^5", -5.19, 54.40)

# leaving out rxn3 
all_rxns = [rxn1, rxn2, rxn4, rxn5, rxn6, rxn7, rxn8, rxn9, rxn10, rxn11, rxn12, rxn13, rxn14,
    rxn15, rxn16, rxn17, rxn18, rxn19, rxn20, rxn21, rxn22, rxn23, rxn24, rxn25, rxn26, rxn27]

@constants R = 8.31446261815324 [unit = u"m_air^3*Pa/K/mol", description = "Universal gas constant"]
@constants PaPerAtm = 101325 [unit = u"Pa/atm", description = "Number of pascals per atmosphere"]
@constants zeroConc = 0.0 [unit = u"mol/m_air^3", description = "Zero concentration"]

@parameters totalK = 1.0e-5 [unit = u"mol/m_air^3", description = "Total concentration of K"]
@parameters totalCa = 1.0e-5 [unit = u"mol/m_air^3", description = "Total concentration of Ca"]
@parameters totalMg = 1.0e-5 [unit = u"mol/m_air^3", description = "Total concentration of Mg"]
@parameters totalNH = 1.0e-5 [unit = u"mol/m_air^3", description = "Total concentration of NH"]
@parameters totalNa = 1.0e-5 [unit = u"mol/m_air^3", description = "Total concentration of Na"]
@parameters totalSO4 = 1.0e-5 [unit = u"mol/m_air^3", description = "Total concentration of SO4"]
@parameters totalNO3 = 1.0e-5 [unit = u"mol/m_air^3", description = "Total concentration of NO3"]
@parameters totalCl = 1.0e-5 [unit = u"mol/m_air^3", description = "Total concentration of Cl"]
@parameters totalH = 1.0e-5 [unit = u"mol/m_air^3", description = "Total concentration of H"]

@variables R_1 = 1.0 [description = "Total sulfate ratio from Section 3.1 of Fountoukis and Nenes (2007)"]
@variables R_2 = 1.0 [description = "Crustal species and sodium ratio from Section 3.1 of Fountoukis and Nenes (2007)"]
@variables R_3 = 1.0 [description = "Crustal species ratio from Section 3.1 of Fountoukis and Nenes (2007)"]

# Simplify activity coefficient for debugging
# logγ₁₂(s::Salt) = 0.5 * I



statevars = [all_solids; all_ions; all_gases; I; W];
params = [T, totalK, totalCa, totalMg, totalNH, totalNa, totalSO4, totalNO3, totalCl, totalH, RH];

function subsystem(rxns, kept_vars)
    eqs = [
        # Add in the reactions.
        [equation(rxn) for rxn ∈ rxns] 
        
        0 ~ SO4_g / SO4_aq # All SO4 immediately goes to aerosol phase as per Section 3.3 (item 1) of Fountoukis and Nenes (2007).
        CaSO4_s ~ 0 # No CaSO4 in solid phase as per Section 3.3 (item 4) of Fountoukis and Nenes (2007).

        # Ratios from Section 3.1
        # R_1 ~ (totalNH + totalCa + totalK + totalMg + totalNa) / totalSO4
        # R_2 ~ (totalCa + totalK + totalMg + totalNa) / totalSO4
        # R_3 ~ (totalCa + totalK + totalMg) / totalSO4

        # Mass balances
        totalK ~ K_aq * W + KHSO4_s + 2K2SO4_s + KNO3_s + KCl_s

        totalCa ~ Ca_aq * W + CaSO4_s + CaNO32_s + CaCl2_s

        totalMg ~ Mg_aq * W + MgSO4_s + MgNO32_s + MgCl2_s

        totalNH ~ NH4_aq * W + NH3_aq * W + NH3_g / (R * T) * PaPerAtm + NH4HSO4_s + 
                    2NH42SO4_s + 3NH43HSO42_s + NH4Cl_s + NH4NO3_s

        totalNa ~ Na_aq * W + NaHSO4_s + 2Na2SO4_s + NaCl_s + NaNO3_s

        totalSO4 ~ (SO4_aq + HSO4_aq) * W + SO4_g / (R * T) * PaPerAtm +
                    KHSO4_s + NaHSO4_s + NH4HSO4_s + CaSO4_s + Na2SO4_s + NH42SO4_s + 
                    2NH43HSO42_s + K2SO4_s + MgSO4_s

        totalNO3 ~ NO3_aq * W + HNO3_aq * W + HNO3_g / (R * T) * PaPerAtm +  NH4NO3_s + NaNO3_s

        totalCl ~ Cl_aq * W + HCl_aq * W + HCl_g / (R * T) * PaPerAtm + 
            NH4Cl_s + NaCl_s + 2CaCl2_s + KCl_s + 2MgCl2_s

        totalH ~ H_aq * W + HNO3_g / (R * T) * PaPerAtm + HCl_g / (R * T) * PaPerAtm

        # Calculate the ionic strength of the multicomponent solution as described by 
        # Fountoukis and Nenes (2007), between equations 8 and 9: ``I = \\frac{1}{2} \\sum_i m_i z_i^2``
        # Force I to always be positive to avoid attempts to take the square root of a negative number.
        I ~ 1 / 2 * sum([ion.m * ion.z^2 for ion in all_Ions])

        # Charge balance
        0 ~ sum([s.cation.m * s.cation.z * s.ν_cation + s.anion.m * s.anion.z * s.ν_anion for s in all_salts])

        # Water content.
        W ~ W_eq16
    ]

    # Set the variables we're not keeping to zero.
    sub_rules = [v => 0.0 for v ∈ setdiff(statevars, kept_vars)]
    eqs = [substitute(eq, Dict(sub_rules)) for eq ∈ eqs]

    # Filter out equations that don't have any of the kept variables.
    function num_variables(eq, kept_vars)
        kept_syms = Symbolics.tosymbol.(kept_vars, (escape=false))
        eq_syms = Symbolics.tosymbol.(get_variables(eq), (escape=false))
        in_vars = eq_syms .∈ Ref(kept_syms)
        sum(in_vars)
    end
    eqs = eqs[[num_variables(eq, kept_vars) > 0 for eq ∈ eqs]]

    # Always work with the absolute value of variables.
    #abs_rules = [v => abs(v) for v ∈ kept_vars]
    #eqs = [substitute(eq, Dict(abs_rules)) for eq ∈ eqs]

    return eqs
end

# Skip unit enforcement for now
# ModelingToolkit.check_units(eqs...) = nothing

regime1_rxns = [rxn5, rxn11, rxn12, rxn13, rxn14, rxn15, rxn16, rxn17, rxn24, rxn26]
# TODO(CT): Missing H2O_g from table 3
regime1_vars = [
    NaHSO4_s, NH4HSO4_s, KHSO4_s, CaSO4_s, 
    Na_aq, NH4_aq, H_aq, HSO4_aq, SO4_aq, NO3_aq, Cl_aq, Ca_aq, K_aq, #H2O_aq, # H2O_g Note: H2O_aq is not a @variable
    NH3_g, NO3_aq, Cl_aq, NH3_aq, HNO3_aq, HCl_aq, # Minor species
    HNO3_g, SO4_g, HCl_g, # Other species not listed in table 3 but still required.
    I, W,
]

regime1eqs = subsystem(regime1_rxns, regime1_vars);
@named regime1sys = NonlinearSystem(regime1eqs, regime1_vars, params)
render(latexify(regime1sys))
guess = ModelingToolkit.get_defaults(regime1sys)
param_defaults = ModelingToolkit.get_defaults(regime1sys)
prob = NonlinearProblem(structural_simplify(regime1sys), guess, param_defaults)
@time sol = solve(prob, NewtonRaphson())

test_rxns = [rxn7, rxn16, rxn17]
test_vars = [
    KCl_s,
    K_aq, Cl_aq, H_aq, OH_aq,
    HCl_g,
    I, W,
]
testeqs = subsystem(test_rxns, test_vars);
testeqs = [substitute(eq, Dict(T => 293.15)) for eq ∈ testeqs]
#testeqs = [substitute(eq, Dict(I => 0.8)) for eq ∈ testeqs]
testeqs = [substitute(eq, Dict(RH => 0.5)) for eq ∈ testeqs]
#testeqs = testeqs[1:end-1]
testeqs = ModelingToolkit.subs_constants(testeqs)
@named testsys = NonlinearSystem(testeqs, test_vars, params)
render(latexify(testsys))

guess = ModelingToolkit.get_defaults(testsys)
param_defaults = ModelingToolkit.get_defaults(testsys)
simplesys = structural_simplify(testsys)
prob = NonlinearProblem(simplesys, guess, param_defaults)
@time sol = solve(prob, NewtonRaphson())

[(u, sol[u]) for u in states(testsys)]
states(simplesys)
render(latexify(simplesys))

aaa = activity(KCl_aqs)
aaa = ModelingToolkit.subs_constants(aaa)
aaa = substitute(aaa, Dict(I => 0.5 * (Cl_aq + H_aq + OH_aq)))
aaa = substitute(aaa, Dict(NO3_aq => 0, NH4_aq => 0, Ca_aq => 0, Na_aq => 0, SO4_aq => 0, H_aq => 0, Mg_aq => 0, OH_aq => 0))
aaa = substitute(aaa, Dict(T => 293.15))

xx = [substitute(aaa, Dict(Cl_aq => cl, K_aq => k)) for cl ∈ 0.01:0.05:5, k ∈ 0.01:0.05:5]
heatmap(0.01:0.05:5, 0.01:0.05:5, Symbolics.value.(xx))

# Tests

# Units from Table 2, last column
@test ModelingToolkit.get_unit(K_eq(rxn5)) == u"mol^3/kg^3"

# Test double activity
@test ModelingToolkit.get_unit(activity([NH3g, HNO3g])) == u"atm^2"

# Check that K_eq decreases as temperature increases.
rxn5_defaults = [rxn5.sys.K⁰ => 24.016, rxn5.sys.C_group => 17.964, rxn5.sys.T₀ => 293.15, rxn5.sys.H_group => -8.423]
@test substitute(K_eq(rxn5), [T => 293.15, rxn5_defaults...]) ≈ 24.016
@test substitute(K_eq(rxn5), [T => 400, rxn5_defaults...]) ≈ 5.545279592576569
@test substitute(K_eq(rxn5), [T => 200, rxn5_defaults...]) ≈ 5429.616126094609

# Test units on all equations
for (i, rxn) in enumerate(all_rxns)
    @test ModelingToolkit.validate(equation(rxn))
end

using ModelingToolkit

@variables a [bounds=(0,1)]
@variables b [bounds=(0,1)]
@variables c [bounds=(0,1)] 

@named sys = NonlinearSystem([
    0.5 ~ a * b,
    1 ~ a + b,
    0.5 ~ a * c,
], [a,b,c], [])

simplesys = structural_simplify(sys)

equations(simplesys)
states(simplesys)

prob = NonlinearProblem(simplesys, [a=>1, b=>2, c=>3], [])
solve(prob, TrustRegion())


@variables t
@variables a(t) b(t) #c(t)
@variables xa(t) xb(t) xtot(t)
@variables Δa(t) Δb(t)
D = Differential(t)

@named sys = ODESystem([
    D(a) ~ 1 - D(b)
    a ~ b
], t, [a,b], [], tspan=(0,1))

@named sys = ODESystem([
    D(a) ~ 1 + Δa * 1e2
    D(b) ~ 0.5 + Δb * 1e2
    xtot ~ a + b
    xtot ~ xa + xb
    2xa ~ xb
    Δa ~ xa - a
    Δb ~ xb - b
    #a ~ xa
    #b ~ xb
], t, [a,b], [], tspan=(0,1))

ssys = structural_simplify(sys)

using DifferentialEquations, Plots
prob = ODEProblem(structural_simplify(sys), [a=>1, b=>1, xa=>1, xb=>1, xtot=>1], (0,1))
sol = solve(prob)
plot(sol)
plot!(sol.t, sol[xb], label="xb(t)")
plot!(sol.t, sol[xtot], label="xtot(t)")

sol[a]
sol[b]
sol[xa]
sol[xb]

using Catalyst
using DifferentialEquations, Plots


rn = @reaction_network begin
    (2*10, 1*10), a <--> b
end

render(latexify(convert(ODESystem, rn)))

prob = ODEProblem(rn, [:a => 1, :b => 0], (0, 2), [])
sol = solve(prob)
plot(sol)