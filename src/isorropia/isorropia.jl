#using Revise
using ModelingToolkit, Catalyst, NonlinearSolve, Unitful, Latexify
using IfElse
using Test
using DifferentialEquations, Plots

# Register the units we'll be using
module MyUnits
using Unitful
@unit m_air "m(air)" MAir 1u"m" false
@unit kg_water "kg(water)" KgWater 1u"kg" false
end
Unitful.register(MyUnits)

# Miscellaneous variables and parameters
@variables t [unit = u"s", description = "Time"]
@variables I(t) = 1.0e-4 [unit = u"mol/kg_water", description = "Ionic strength"]
@variables W(t) = 1.0e-10 [unit = u"kg_water/m_air^3", description = "Aerosol water content"]
for v ∈ (:I, :W)
    eval(:($v = ParentScope($v))) # Keep these as global MTK variables.
end

@parameters T = 293.15 [unit = u"K", description = "Temperature"]
@parameters RH = 0.3 [description = "Relative humidity (expressed on a scale from 0 to 1)"] # unitless
for p ∈ (:T, :RH)
   eval(:($p = ParentScope($p)))
end

@constants tiny_conc=1e-20 [unit = u"mol/m_air^3", description = "Tiny concentration to avoid division by zero"]

"""
A species represents a chemical species in the system.

Chemical species should have an `γ` method which returns the activity coefficient
of the species as specified in Section 2.2 of Fountoukis and Nenes (2007). 
They should also have a `terms` method which returns the variable(s) and
stoichiometry coefficient(s) associated with the species, and also a 
`min_conc` method which returns the minimum concentration among the terms that 
make up the species.
"""
abstract type Species end
activity(s) = reduce(*, [m^ν for (m, ν) ∈ zip(terms(s)...)]) * γ(s)
γ(s::Species) = error("activity coefficient γ not defined for $(typeof(s))")
terms(s::Species) = error("terms not defined for $(typeof(s))")

""" 
The activity coefficient of multiple species is the product of their Activity
coefficients as shown in Table 2 of Fountoukis and Nenes (2007).
"""
γ(s::AbstractVector) = reduce(*, γ.(s))
function terms(s::AbstractVector) 
    tt = terms.(s)
    vcat([t[1] for t ∈ tt]...), vcat([t[2] for t ∈ tt]...)
end
min_conc(s::AbstractArray) = reduce(min, min_conc.(s))

@constants R = 8.31446261815324 [unit = u"m_air^3*Pa/K/mol", description = "Universal gas constant"]
@constants PaPerAtm = 101325 [unit = u"Pa/atm", description = "Number of pascals per atmosphere"]
@constants W_one = 1 [unit = u"kg_water/m_air^3", description = "Unit aerosol water content"]
mol2atm(p) = p / PaPerAtm * R * T

# Load the other files
include("isorropia_aqueous.jl")
include("isorropia_deliquescence.jl")
include("isorropia_solid.jl")
include("isorropia_gas.jl")
mw = Dict(Na_aq => 22.989769, SO4_aq => 96.0636, SO4_g => 96.0636, NH3_aq => 17.03052, NH3_g => 17.03052, NO3_aq => 62.0049, Cl_aq => 35.453, 
        Ca_aq => 40.078, K_aq => 39.0983, Mg_aq => 24.305, H_aq => 1.00784, NH4_aq => 18.04,
        K2SO4_s => 174.259, KNO3_s => 101.1032, CaNO32_s => 164.1, HNO3_g => 63.01, HNO3_aq => 63.01) # g/mol
include("isorropia_water.jl")

@constants tinyrate = 1.e-99 [unit = u"mol/m_air^3/s", description = "Tiny reaction rate constant negative concentrations"]

"""
Define a reaction based on information in Table 2 of Fountoukis and Nenes (2007).
"""
struct Rxn
    reactant
    product
    sys

    function Rxn(reactant, product, K⁰::Number, K⁰units, hgroup::Number, cgroup::Number; name="rxn")
        γr = γ(reactant)
        γp = γ(product)
        ar = activity(reactant)
        ap = activity(product)
        # These are the variables from Fountoukis and Nenes (2007) Table 2
        @constants K⁰ = K⁰ [unit = K⁰units, description = "Equilibrium constant at 298.15 K"] 
        @constants H_group = hgroup [description = "ΔH⁰ / (R * T₀) (unitless)"]
        @constants C_group = cgroup [description = "ΔC⁰ₚ / R (unitless)"]
        @variables K_eq(t) [unit = K⁰units, description = "Equilibrium constant"]
        # These are the transformed variables to turn this into a kinetic reaction rather than an equilibrium reaction.
        # To do so we need to choose a base reaction rate constant. (We choose it to be 1.0e-15 here.)
        @constants fakerate = 1e6 [unit = u"mol/m_air^3/s", description = "Fake reaction rate constant to convert equilibrium reaction into kinetic reaction"]
        @constants γrkludge = 1 [unit = ModelingToolkit.get_unit(1/ar), description = "Unit conversion for γ_r"]
        @constants γpkludge = 1 [unit = ModelingToolkit.get_unit(1/ap), description = "Unit conversion for γ_p"]
        @constants Keqkludge = 1 [unit = ModelingToolkit.get_unit(ar/ap), description = "Unit conversion for K_eq"]
        @variables k_fwd(t)=1 [unit = ModelingToolkit.get_unit(γr/ar*fakerate), description = "Forward reaction rate constant"]
        @variables k_rev(t)=1 [unit = ModelingToolkit.get_unit(γp/ap*fakerate), description = "Reverse reaction rate constant"]
        @constants tinyk_fwd=1e-20 [unit = ModelingToolkit.get_unit(γr/ar*fakerate), description = "Minimum forward reaction rate constant"]
        @constants tinyk_rev=1e-20 [unit = ModelingToolkit.get_unit(γp/ap*fakerate), description = "Minimum reverse reaction rate constant"]
        @constants unit_kfwd=1 [unit = ModelingToolkit.get_unit(γr/ar*fakerate), description = "Minimum forward reaction rate constant"]
        @constants unit_krev=1 [unit = ModelingToolkit.get_unit(γp/ap*fakerate), description = "Unit conversion for k_rev"]
        @variables eq_ratio(t)=1 [unit = ModelingToolkit.get_unit(k_fwd/k_rev), description = "Unitless equilibrium ratio of product to reactant"]
        @variables rate_fwd = 1 [unit = u"mol/m_air^3/s", description = "Forward reaction rate law"]
        @variables rate_rev = 1 [unit = u"mol/m_air^3/s", description = "Reverse reaction rate law"]
        @variables γ_r(t) [unit=ModelingToolkit.get_unit(γr), description="Activity coefficient of reactant"]
        @variables γ_p(t) [unit=ModelingToolkit.get_unit(γp), description="Activity coefficient of product"]
        @constants min_γr = 1e-5 [unit=ModelingToolkit.get_unit(γr), description="Minimum activity coefficient of reactant"]
        @constants min_γp = 1e-5 [unit=ModelingToolkit.get_unit(γp), description="Minimum activity coefficient of product"]   
        if (typeof(product) <: SaltLike) && (reactant isa Solid) # Deliquescence
            krev = ifelse(min_conc(product) > tiny_conc, 
                f_deliquescence[product] * γ_p* γpkludge * fakerate, 
                #γ_p * γpkludge * fakerate, 
                #tinyk_rev, 
                tinyk_rev)
        else
            krev = ifelse(min_conc(product) > tiny_conc, γ_p* γpkludge * fakerate, tinyk_rev)
        end
        kfwd = ifelse(min_conc(reactant) > tiny_conc, K_eq * γ_r * Keqkludge * γrkludge * fakerate, tinyk_fwd)
        sys = NonlinearSystem([
            K_eq ~ K⁰ * exp(-H_group * (T₀ / T - 1) - C_group * (1 + log(T₀ / T) - T₀ / T))
            γ_r ~ max(min_γr, γ(reactant)) 
            γ_p ~ max(min_γp, γ(product)) 
            k_fwd ~ max(min(kfwd / krev * unit_krev, unit_kfwd*1e10), unit_kfwd*1e-10)
            k_rev ~ unit_krev
            # Only react if there is something to react.
            #rate_fwd ~ ifelse(min_conc(reactant) > tiny_conc, K_eq * Keqkludge * γrkludge * fakerate * ar, tinyrate)
            #rate_rev ~ rrev
            #eq_ratio ~ k_fwd / k_rev
        ], [K_eq, γ_r, γ_p, k_fwd, k_rev, rate_fwd, rate_rev, eq_ratio], [T₀, K⁰, H_group, C_group]; name=name)
        new(reactant, product, sys)
    end
end

# Equation 5: Equilibrium constant
@constants T₀ = 293.15 [unit = u"K", description = "Standard temperature"]

"""
Assemble an equation for a reaction based on Table 2 of Fountoukis and Nenes (2007), where
the left-hand side is the equilibrium constant and the right-hand side is the activity.
"""
function reactions(r::Rxn)
    rterms = terms(r.reactant)
    pterms = terms(r.product)
    fwd = Reaction(r.sys.k_fwd,  rterms[1], pterms[1], rterms[2], pterms[2])#, only_use_rate=true)
    rev = Reaction(r.sys.k_rev, pterms[1], rterms[1], pterms[2], rterms[2])#, only_use_rate=true)
    return [fwd, rev]
end


# NOTE: Assuming that H_group and C_group are zero when they are left out of Table 2. 
@named rxn1 = Rxn(CaNO32s, CaNO32_aqs, 6.067e5, u"mol^3/kg_water^3", -11.299, 0.0)
@named rxn2 = Rxn(CaCl2s, CaCl2_aqs, 7.974e11, u"mol^3/kg_water^3", -14.087, 0.0)
@named rxn3 = Rxn(CaSO4s, CaSO4_aqs, 4.319e-5, u"mol^2/kg_water^2", 0.0, 0.0)
@named rxn4 = Rxn(K2SO4s, K2SO4_aqs, 1.569e-2, u"mol^3/kg_water^3", -9.589, 45.807)
@named rxn5 = Rxn(KHSO4s, KHSO4_aqs, 24.016, u"mol^2/kg_water^2", -8.423, 17.964)
@named rxn6 = Rxn(KNO3s, KNO3_aqs, 0.872, u"mol^2/kg_water^2", 14.075, 19.388)
@named rxn7 = Rxn(KCls, KCl_aqs, 8.680, u"mol^2/kg_water^2", -6.167, 19.953)
@named rxn8 = Rxn(MgSO4s, MgSO4_aqs, 1.079e5, u"mol^2/kg_water^2", 36.798, 0.0)
@named rxn9 = Rxn(MgNO32s, MgNO32_aqs, 2.507e15, u"mol^3/kg_water^3", -8.754, 0.0)
@named rxn10 = Rxn(MgCl2s, MgCl2_aqs, 9.557e21, u"mol^3/kg_water^3", -1.347, 0.0)
@named rxn11 = Rxn(HSO4_ion, [H_ion, SO4_ion], 1.015e-2, u"mol/kg_water", 8.85, 25.14)
@named rxn12 = Rxn(NH3g, NH3_ion, 5.764e1, u"mol/kg_water/atm", 13.79, -5.39)
@named rxn13 = Rxn([NH3_ion, H2Oaq], [NH4_ion, OH_ion], 1.805e-5, u"mol/kg_water", -1.50, 26.92)
@named rxn14 = Rxn(HNO3g, HNO3_aqs, 2.511e6, u"mol^2/kg_water^2/atm", 29.17, 16.83)
@named rxn15 = Rxn(HNO3g, HNO3_ion, 2.1e5, u"mol/kg_water/atm", 29.17, 16.83)
@named rxn16 = Rxn(HClg, [H_ion, Cl_ion], 1.971e6, u"mol^2/kg_water^2/atm", 30.20, 19.91)
@named rxn17 = Rxn(HClg, HCl_ion, 2.5e3, u"mol/kg_water/atm", 30.20, 19.91)
@named rxn18 = Rxn(H2Oaq, [H_ion, OH_ion], 1.010e-14, u"mol^2/kg_water^2", -22.52, 26.92)
@named rxn19 = Rxn(Na2SO4s, Na2SO4_aqs, 4.799e-1, u"mol^3/kg_water^3", 0.98, 39.75)
@named rxn20 = Rxn(NH42SO4s, NH42SO4_aqs, 1.87e0, u"mol^3/kg_water^3", -2.65, 38.57)
@named rxn21 = Rxn(NH4Cls, [NH3g, HClg], 1.086e-16, u"atm^2", -71.00, 2.40)
@named rxn22 = Rxn(NaNO3s, NaNO3_aqs, 1.197e1, u"mol^2/kg_water^2", -8.22, 16.01)
@named rxn23 = Rxn(NaCls, NaCl_aqs, 3.766e1, u"mol^2/kg_water^2", -1.56, 16.90)
@named rxn24 = Rxn(NaHSO4s, NaHSO4_aqs, 2.413e4, u"mol^2/kg_water^2", 0.79, 14.75)
@named rxn25 = Rxn(NH4NO3s, [NH3g, HNO3g], 4.199e-17, u"atm^2", -74.375, 6.025)
@named rxn26 = Rxn(NH4HSO4s, NH4HSO4_aqs, 1.383e0, u"mol^2/kg_water^2", -2.87, 15.83)
@named rxn27 = Rxn(NH43HSO42s, NH43HSO42_aqs, 2.972e1, u"mol^5/kg_water^5", -5.19, 54.40)

all_rxns = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6, rxn7, rxn8, rxn9, rxn10, rxn11, rxn12, rxn13, rxn14,
    rxn15, rxn16, rxn17, rxn18, rxn19, rxn20, rxn21, rxn22, rxn23, rxn24, rxn25, rxn26, rxn27]

@test ModelingToolkit.get_unit(mol2atm(SO4_g)) == u"atm"

statevars = [all_solids; all_ions; all_gases; I; W]
ps = [T, RH, H2O_aq, metastable]

@named IW = NonlinearSystem([     
        # Calculate the ionic strength of the multicomponent solution as described by 
        # Fountoukis and Nenes (2007), between equations 8 and 9: ``I = \\frac{1}{2} \\sum_i m_i z_i^2``
        # Force I to always be positive to avoid attempts to take the square root of a negative number.
        I ~ max(1.0e-20*I_one, 1 / 2 * sum([ion.m * ion.z^2 for ion in all_Ions] / W))

        # Water content.
        W ~ max(1.0e-10*W_one, W_eq16)
], [I; W], [])

x = vcat([reactions(x) for x in all_rxns]...)
@parameters so4rate=100 [unit = u"s^-1", description = "Rate of SO4 to aerosol phase (pseudo-instantaneous)"]
push!(x, Reaction(ifelse(SO4_g > tiny_conc, so4rate, tinyrate/unit_conc), [SO4_g], [SO4_aq])) # All SO4 immediately goes to aerosol phase as per Section 3.3 (item 1) of Fountoukis and Nenes (2007).
xx = ReactionSystem(x, t, statevars, [so4rate; ps]; 
    systems=[IW, DRH, [rxn.sys for rxn ∈ all_rxns]...], 
    checks=false, combinatoric_ratelaws=false, name=:xx)
#render(latexify(xx))
testsys = convert(ODESystem, xx)
#render(latexify(testsys))

@test all([ModelingToolkit.check_units(eq) for eq in testsys.eqs])

pp = structural_simplify(testsys)
u₀ = ModelingToolkit.get_defaults(pp)
u₀[RH] = 0.1
u₀[HNO3_g] = 2 / 1e6 / mw[HNO3_g]
#u₀[CaNO32_s] = 2 / 1e6 / mw[CaNO32_s]
prob = ODEProblem(pp, u₀, (0.0, 10000.0), u₀)
@time sol = solve(prob, Rosenbrock23())

function plotvars(sol, vars; kwargs...)
    p1 = plot(; kwargs...)
    for v ∈ vars
        plot!(p1, sol[t], sol[v], label=v)
    end
    p1
end

plot([
plotvars(sol, [I]; title="I")
plotvars(sol, [W]; title="W")
plotvars(sol, states(testsys)[[occursin("aq", string(v)) && (string(v) != "H2O_aq(t)") for v ∈ states(testsys)]];
    title="Aqueous species")
plotvars(sol, states(testsys)[[occursin("_s", string(v)) for v ∈ states(testsys)]]; title="Solids")
plotvars(sol, states(testsys)[[occursin("_g", string(v)) for v ∈ states(testsys)]]; title="Gases")
plotvars(sol, states(testsys)[[occursin("k_fwd", string(v)) for v ∈ states(testsys)]]; title="Forward rates", yscale=:log10)
plotvars(sol, states(testsys)[[occursin("k_rev", string(v)) for v ∈ states(testsys)]]; title="Reverse rates", yscale=:log10)
]..., size=(1500, 1000))


plot(
    plot(sol[t], sol[DRH.f_CaNO32], 
        ylabel="f_CaNO32", xlabel="time (s)", label=:none),
    plot(sol[t], sol[CaNO32_s] * 1e6 * mw[CaNO32_s], 
        ylabel="CaNO32_s", xlabel="time (s)", label=:none),
    plot(sol[t], sol[NO3_aq] * 1e6 * mw[NO3_aq], 
        ylabel="NO3_aq", xlabel="time (s)", label=:none),
    plot(sol[t], sol[HNO3_aq] * 1e6 * mw[HNO3_aq], 
        ylabel="HNO3_aq", xlabel="time (s)", label=:none),
    plot(sol[t], sol[HNO3_g] * 1e6 * mw[HNO3_g], 
        ylabel="HNO3_g", xlabel="time (s)", label=:none),
    plot(sol[t], sol[Ca_aq] * 1e6 * mw[Ca_aq], 
        ylabel="Ca_aq", xlabel="time (s)", label=:none),
    plot(sol[t], sol[rxn1.sys.k_fwd], 
        ylabel="k_fwd", xlabel="time (s)", label=:none, yscale=:log10),
    plot(sol[t], sol[rxn1.sys.k_rev], 
        ylabel="k_rev", xlabel="time (s)", label=:none, yscale=:log10),
)

plot([plot(sol[t], sol[ion], ylabel=ion, xlabel="time", label=:none) for ion in all_ions]..., size=(1000, 800))
plot([plot(sol[t], sol[gas], ylabel=gas, xlabel="time", label=:none) for gas in all_gases]..., size=(1000, 800))
plot([plot(sol[t], sol[solid], ylabel=solid, xlabel="time", label=:none) for solid in all_solids]..., size=(1000, 800))


@testset "Mass balances" begin 
    totalK(sol) = [s * sol[x] for (s, x) in [(1, K_aq), (1, KHSO4_s), (2, K2SO4_s), (1, KNO3_s), (1, KCl_s)]]
    totalCa(sol) = [s * sol[x] for (s, x) in [(1, Ca_aq), (1, CaSO4_s),  (1, CaNO32_s), (1, CaCl2_s)]]
    totalMg(sol) = [s * sol[x] for (s, x) in [(1, Mg_aq), (1, MgSO4_s), (1, MgNO32_s), (1, MgCl2_s)]]
    totalNH(sol) = [s * sol[x] for (s, x) in [(1, NH4_aq), (1, NH3_aq), (1, NH3_g), (1, NH4HSO4_s), 
                (2, NH42SO4_s), (3, NH43HSO42_s), (1, NH4Cl_s), (1, NH4NO3_s)]]
    totalNa(sol) = [s * sol[x] for (s, x) in [(1, Na_aq), (1, NaHSO4_s), (2, Na2SO4_s), (1, NaCl_s), (1, NaNO3_s)]]
    totalSO4(sol) = [s * sol[x] for (s, x) in [(1, SO4_aq), (1, HSO4_aq), (1, SO4_g),
                    (1, KHSO4_s), (1, NaHSO4_s), (1, NH4HSO4_s), (1, CaSO4_s), (1, Na2SO4_s), (1, NH42SO4_s), 
                    (2, NH43HSO42_s), (1, K2SO4_s), (1, MgSO4_s)]]
    totalNO3(sol) = [s * sol[x] for (s, x) in [(1, NO3_aq), (1, HNO3_aq), (1, HNO3_g),  (1, NH4NO3_s), (1, NaNO3_s)]]
    totalCl(sol) = [s * sol[x] for (s, x) in [(1, Cl_aq), (1, HCl_aq), (1, HCl_g), (1, NH4Cl_s),
                        (1, NaCl_s), (2, CaCl2_s), (1, KCl_s), (2, MgCl2_s)]]
    totalH(sol) = [s * sol[x] for (s, x) in [(1, H_aq), (1, HNO3_g), (1, HCl_g)]]

    for tot ∈ [totalK, totalCa, totalMg, totalNH, totalNa, totalSO4, totalNO3, totalCl, totalH]
        @testset "$(string(tot))" begin
            total = sum(hcat(tot(sol)...), dims=2)
            @test (maximum(abs.(total .- sum(total)/length(total)))) < 1.e-10
        end
    end
end

prob = SteadyStateProblem(pp, u₀, u₀)
@time sol = solve(prob, DynamicSS(Rosenbrock23()))
sol[CaNO32_s]  * 1e6 * mw[CaNO32_s]

# Fountoukis and Nenes (2007) Figure 6
u₀ = ModelingToolkit.get_defaults(pp)
ics = Dict([Na_aq => 0, SO4_g => 10, NH3_g => 3.4, HNO3_g => 2, Cl_aq => 0, 
    Ca_aq => 0.4, K_aq => 0.33, Mg_aq => 1e-10]) # ug/m3
for k ∈ keys(ics)
    u₀[k] = ics[k] / 1e6 / mw[k] # ug/m3 / (1e6 ug/g) / g/mol = mol/m3
end
u₀[H_aq] = 2*u₀[SO4_aq] + u₀[NO3_aq] + u₀[Cl_aq]
for g ∈ setdiff(all_gases, ics)
    u₀[g] = 1e-10
end
u₀[metastable] = 0

RHs = [10, 25, 40, 55, 65, 70, 75, 80, 85, 90] ./ 100.0
sols = []
for rh in RHs
    u₀[RH] = rh
    #prob = ODEProblem(pp, u₀, (0.0, 1.0), u₀)
    #sol = solve(prob, Rosenbrock23())
    prob = SteadyStateProblem(pp, u₀, u₀)
    sol = solve(prob, DynamicSS(Rosenbrock23()))
    push!(sols, sol)
end

plot(
    plot(RHs, [sols[i][W] * 1e9 for i ∈ 1:length(RHs)], ylim=(0, 50),
        ylabel="H2O (ug/m3)", xlabel="Relative humidity (%)", label=:none),
    plot(RHs, [sols[i][K_aq] * 1e6 * mw[K_aq] for i ∈ 1:length(RHs)], ylim=(0,1.2), 
        ylabel="K+(aq) (ug/m3)", xlabel="Relative humidity (%)", label=:none),
    plot(RHs, [sols[i][NH4_aq] * 1e6 * mw[NH4_aq] for i ∈ 1:length(RHs)], ylim=(0,3.0),
        ylabel="NH4(aq) (ug/m3)", xlabel="Relative humidity (%)", label=:none),
    plot(RHs, [sols[i][NO3_aq] * 1e6 * mw[NO3_aq] for i ∈ 1:length(RHs)], #ylim=(0, 0.25),
        ylabel="NO3(aq) (ug/m3)", xlabel="Relative humidity (%)", label=:none),
)

plot(
plot(RHs, [sols[i][DRH.f_K2SO4] for i ∈ 1:length(RHs)], 
    ylabel="f_K2SO4", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][K2SO4_s] * 1e6 * mw[K2SO4_s] for i ∈ 1:length(RHs)], 
    ylabel="K2SO4_s", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][SO4_aq] * 1e6 * mw[SO4_aq] for i ∈ 1:length(RHs)], 
    ylabel="SO4_aq", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][K_aq] * 1e6 * mw[K_aq] for i ∈ 1:length(RHs)], 
    ylabel="K_aq", xlabel="Relative humidity (%)", label=:none),
    plot(RHs, [sols[i][rxn4.sys.k_fwd] for i ∈ 1:length(RHs)], 
    ylabel="k_fwd", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][rxn4.sys.k_rev] for i ∈ 1:length(RHs)], 
    ylabel="k_rev", xlabel="Relative humidity (%)", label=:none),
)

plot(
plot(RHs, [sols[i][DRH.f_CaNO32] for i ∈ 1:length(RHs)], 
    ylabel="f_CaNO32", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][CaNO32_s] * 1e6 * mw[CaNO32_s] for i ∈ 1:length(RHs)], 
    ylabel="CaNO32_s", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][NO3_aq] * 1e6 * mw[NO3_aq] for i ∈ 1:length(RHs)], 
    ylabel="NO3_aq", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][Ca_aq] * 1e6 * mw[Ca_aq] for i ∈ 1:length(RHs)], 
    ylabel="Ca_aq", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][rxn1.sys.k_fwd] for i ∈ 1:length(RHs)], 
    ylabel="k_fwd", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][rxn1.sys.k_rev] for i ∈ 1:length(RHs)], 
    ylabel="k_rev", xlabel="Relative humidity (%)", label=:none),
)

plot(
plot(RHs, [sols[i][DRH.f_KNO3] for i ∈ 1:length(RHs)], 
    ylabel="f_KNO3", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][KNO3_s] * 1e6 * mw[KNO3_s] for i ∈ 1:length(RHs)], 
    ylabel="KNO3_s", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][NO3_aq] * 1e6 * mw[NO3_aq] for i ∈ 1:length(RHs)], 
    ylabel="NO3_aq", xlabel="Relative humidity (%)", label=:none),
plot(RHs, [sols[i][K_aq] * 1e6 * mw[K_aq] for i ∈ 1:length(RHs)], 
    ylabel="K_aq", xlabel="Relative humidity (%)", label=:none),
    plot(RHs, [sols[i][rxn6.sys.k_fwd] for i ∈ 1:length(RHs)], 
    ylabel="k_fwd", xlabel="Relative humidity (%)", label=:none, yscale=:log10),
plot(RHs, [sols[i][rxn6.sys.k_rev] for i ∈ 1:length(RHs)], 
    ylabel="k_rev", xlabel="Relative humidity (%)", label=:none, yscale=:log10),
)


plot([plot(RHs, [sols[i][ion] for i ∈ 1:length(RHs)], 
    ylabel=ion, xlabel="Relative humidity (%)", label=:none) for ion in all_ions]..., size=(1000, 800))

plot([plot(RHs, [sols[i][gas] for i ∈ 1:length(RHs)], 
    ylabel=gas, xlabel="Relative humidity (%)", label=:none) for gas in all_gases]..., size=(1000, 800))

plot([plot(RHs, [sols[i][solid] for i ∈ 1:length(RHs)], 
    ylabel=solid, xlabel="Relative humidity (%)", label=:none) for solid in all_solids]..., size=(1000, 800))




uu = Dict([(k, 0.01) for k ∈ all_ions])
x = substitute(W_eq16, uu)
rhs = 0.01:0.01:1.0
plot(rhs, [Symbolics.value(substitute(x, Dict(RH => w, unit_molality=>1))) for w ∈ rhs])


show([(x => sol[x][1]) for x in states(testsys)])
[(x => eval(:(sol[$x.sys.eq_ratio][1]))) for x in Symbol.("rxn".*string.(1:27))]



syms = [Symbol("rxn$i") for i in 1:27]
vs = [eval(:(sol[$(s).sys.eq_ratio])) for s ∈ syms]
show(vcat(ModelingToolkit.value.(vs)...))

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
@named regime1rsys = ReactionSystem(regime1eqs, t, regime1_vars, params)
#@named regime1sys = NonlinearSystem(regime1eqs, regime1_vars, params)
render(latexify(regime1rsys))
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