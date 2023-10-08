using EarthSciMLBase
using ModelingToolkit, Catalyst, DifferentialEquations, Unitful
using Test
using Plots

include(joinpath(@__DIR__, "../../src/isorropia/isorropia.jl"))

@variables t [unit = u"s", description = "Time"]

model = Isorropia(t);

sys = structural_simplify(get_mtk(model))

@test all([ModelingToolkit.check_units(eq) for eq in equations(get_mtk(model))])


mw = Dict(Na_aq => 22.989769, SO4_aq => 96.0636, SO4_g => 96.0636, NH3_aq => 17.03052, NH3_g => 17.03052,
    NO3_aq => 62.0049, Cl_aq => 35.453,
    Ca_aq => 40.078, K_aq => 39.0983, Mg_aq => 24.305, H_aq => 1.00784, NH4_aq => 18.04, HCl_g => 36.46,
    K2SO4_s => 174.259, KNO3_s => 101.1032, CaNO32_s => 164.1, HNO3_g => 63.01, HNO3_aq => 63.01,
    KHSO4_s => 136.169, KCl_s => 74.5513, NH4NO3_s => 80.043) # g/mol

# Molecules for checking mass balance
K_molecs = [(1, K_aq), (1, KHSO4_s), (2, K2SO4_s), (1, KNO3_s), (1, KCl_s)]
Ca_molecs = [(1, Ca_aq), (1, CaSO4_s), (1, CaNO32_s), (1, CaCl2_s)]
Mg_molecs = [(1, Mg_aq), (1, MgSO4_s), (1, MgNO32_s), (1, MgCl2_s)]
NH_molecs = [(1, NH4_aq), (1, NH3_aq), (1, NH3_g), (1, NH4HSO4_s),
    (2, NH42SO4_s), (3, NH43HSO42_s), (1, NH4Cl_s), (1, NH4NO3_s)]
Na_molecs = [(1, Na_aq), (1, NaHSO4_s), (2, Na2SO4_s), (1, NaCl_s), (1, NaNO3_s)]
SO4_molecs = [(1, SO4_aq), (1, HSO4_aq), (1, SO4_g),
    (1, KHSO4_s), (1, NaHSO4_s), (1, NH4HSO4_s), (1, CaSO4_s), (1, Na2SO4_s), (1, NH42SO4_s),
    (2, NH43HSO42_s), (1, K2SO4_s), (1, MgSO4_s)]
NO3_molecs = [(1, NO3_aq), (1, HNO3_aq), (1, HNO3_g), (1, NH4NO3_s), (1, NaNO3_s)]
Cl_molecs = [(1, Cl_aq), (1, HCl_aq), (1, HCl_g), (1, NH4Cl_s),
    (1, NaCl_s), (2, CaCl2_s), (1, KCl_s), (2, MgCl2_s)]
H_molecs = [(1, H_aq), (1, HNO3_g), (1, HCl_g)]


defaults = ModelingToolkit.get_defaults(sys)
u₀ = Dict{Any,Float64}([s => 1.0e-15 for s ∈ states(sys)])
u₀[CaNO32_s] = 0.8 / 1e6 / mw[CaNO32_s]
u₀[Ca_aq] = 0.8 / 1e6 / mw[Ca_aq]
u₀[NO3_aq] = 10.0 / 1e6 / mw[NO3_aq]

p = Dict{Any,Float64}([p => defaults[p] for p ∈ parameters(sys)])
p[RH] = 0.30
prob = ODEProblem(sys, u₀, (0.0, 100.0), p)
# Need low tolerance for mass balance checks to pass.
@time sol = solve(prob, Rosenbrock23(), abstol=1e-12, reltol=1e-12)

eq = []
for i in eachindex(sol.u)
    u = [x => sol[x][i] for x in [I, W, states(sys)...]]
    push!(u, T => sol[T])
    γ_p = substitute(ModelingToolkit.subs_constants(γ(CaNO32_aqs)), u)
    a_r = substitute(ModelingToolkit.subs_constants(activity(CaNO32s)), u)
    a_p = substitute(ModelingToolkit.subs_constants(activity(CaNO32_aqs)), u)
    keq = sol[rxn1.sys.K_eq][i]
    @info :a_p => a_p, :γ_p => γ_p, :Ca_aq => sol[Ca_aq][i], :a_r => a_r, :keq => keq, :ratio => a_p / a_r / keq
    push!(eq, a_p / a_r / keq)
end

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
    plot(sol[t], sol[rxn1.sys.k_rev],
        ylabel="k_rev", xlabel="time (s)", label=:none),
    plot(sol[t], sol[W] * 1e9,
        ylabel="W (ug/m3)", xlabel="time (s)", label=:none),
)

plot([plot(sol[t], sol[ion], ylabel=ion, xlabel="time", label=:none) for ion in all_ions]..., size=(1000, 800))
plot([plot(sol[t], sol[gas], ylabel=gas, xlabel="time", label=:none) for gas in all_gases]..., size=(1000, 800))
plot([plot(sol[t], sol[solid], ylabel=solid, xlabel="time", label=:none) for solid in all_solids]..., size=(1000, 800))

@testset "Mass balances" begin
    names = [:K, :Ca, :Mg, :NH, :Na, :SO4, :NO3, :Cl, :H]
    for (i, molecs) ∈ enumerate([K_molecs, Ca_molecs, Mg_molecs, NH_molecs, Na_molecs, SO4_molecs, NO3_molecs, Cl_molecs, H_molecs])
        tot(sol) = [s * sol[x] for (s, x) in molecs]
        @testset "$(names[i])" begin
            total = sum(hcat(tot(sol)...), dims=2)
            if names[i] ∈ [:NO3, :H]
                @test_broken (maximum(abs.(total .- sum(total) / length(total)))) < 1.e-10
            else
                @test (maximum(abs.(total .- sum(total) / length(total)))) < 1.e-10
            end
        end
    end
end

# TOOD(CT): ODEProblem and SteadyStateProblem don't give the same result.
prob = SteadyStateProblem(sys, u₀, p)
@time sol_ss = solve(prob, DynamicSS(Rosenbrock23()), abstol=1e-12, reltol=1e-12)

plot(
    bar(max.((1.e-22), sol.u[end]), yscale=:log10, legend=:none, title="Not steady state", ylim=(1e-10, 1e-6)),
    bar(max.((1.e-22), sol_ss.u), yscale=:log10, legend=:none, title="Steady state", ylim=(1e-10, 1e-6))
)


# Check equilibrium for the first equation
defaults = ModelingToolkit.get_defaults(sys)
u₀ = Dict{Any,Float64}([s => 1.0e-20 for s ∈ states(sys)])
u₀[CaNO32_s] = 0.8 / 1e6 / mw[CaNO32_s]
u₀[Ca_aq] = 0.8 / 1e6 / mw[Ca_aq]
u₀[NO3_aq] = 1.6 / 1e6 / mw[Ca_aq]

p = Dict{Any,Float64}([p => defaults[p] for p ∈ parameters(sys)])
p[RH] = 0.02

# prob = SteadyStateProblem(sys, u₀, p)
# @time sol = solve(prob, DynamicSS(Rosenbrock23()), abstol=1e-12, reltol=1e-12)

prob = ODEProblem(sys, u₀, (0.0, 30.0), p)
@time sol = solve(prob, Rosenbrock23(), abstol=1e-12)

u = [x => sol[x][end] for x in [I, W, T, states(sys)...]]
ca2plus = sol[Ca_aq][end] / sol[W][end] # mol/kg_water
no3minus = sol[NO3_aq][end] / sol[W][end] # mol/kg_water

γaq = substitute(ModelingToolkit.subs_constants(
        exp(logγ₁₂(CaNO32_aqs))^(CaNO32_aqs.ν_cation + CaNO32_aqs.ν_anion)), u)

# First activity calculated using formulas from paper.
aq_activity1 = ca2plus * no3minus^2 * γaq

# Second activity calculated using functions with water conversions.
aq_activity2 = substitute(ModelingToolkit.subs_constants(activity(CaNO32_aqs)), u)

@test ModelingToolkit.value(aq_activity1) ≈ ModelingToolkit.value(aq_activity2)

keq = rxn1.sys.K⁰ * exp(-rxn1.sys.H_group * (T₀ / T - 1) -
                        rxn1.sys.C_group * (1 + log(T₀ / T) - T₀ / T))
eq_const = substitute(ModelingToolkit.subs_constants(keq), u)

solid_activity = substitute(ModelingToolkit.subs_constants(activity(CaNO32s)), u)
@test solid_activity ≈ 1.0

# Theoretically the ratio of the activities should be equal to the equilibrium constant
# at the end of the simulation, but it seems like in practice sometimes the concentration on one
# side or the other gets too small and the equilibrium constant blows up.
@test isapprox(ModelingToolkit.value(aq_activity2) /
      ModelingToolkit.value(solid_activity) / ModelingToolkit.value(eq_const), 1.0, atol=1e-2)

# Derivative of activity with respect to concentration should be positive
da_daq = substitute(ModelingToolkit.subs_constants(
        expand_derivatives(Differential(Ca_aq)(activity(CaNO32_aqs)))), [I_one => 1.0, T₀₂ => 273.15, c_1 => 0.005, u...])
@test ModelingToolkit.value(da_daq) > 0.0

# The derivative of our equilibrium ratio (the ratio of our equilibrium expression to one)
# is positive for the aqueous concentration and zero for the solid concentration, because
# the rate of the forward reaction (solid to aqueous) doesn't depend on the concentration of the solid,
# but rate of the reverse reaction (aqueous to solid) does depend on concentration of the aqueous salt.
k_expr = ModelingToolkit.subs_constants(activity(CaNO32_aqs) / activity(CaNO32s) / keq)
dk_daq = expand_derivatives(Differential(Ca_aq)(k_expr))

@test substitute(dk_daq, u) > 0

dk_ds = expand_derivatives(Differential(CaNO32_s)(k_expr))
@test substitute(dk_ds, u₀) == 0

@test ModelingToolkit.get_unit(activity(CaNO32_aqs) / activity(CaNO32s) / keq) isa Unitful.FreeUnits{(),NoDims,nothing}







# Fountoukis and Nenes (2007) Figure 6
defaults = ModelingToolkit.get_defaults(sys)
ics = Dict([Na_aq => 0, SO4_g => 10, NH3_g => 3.4, HNO3_g => 2, HCl_g => 0,
   Ca_aq => 0.4, K_aq => 0.33, Mg_aq => 1e-20]) # ug/m3
for k ∈ keys(ics)
   u₀[k] = ics[k] / 1e6 / mw[k] # ug/m3 / (1e6 ug/g) / g/mol = mol/m3
end
u₀[H_aq] = 2 * u₀[SO4_aq]

p = Dict{Any,Float64}([p => defaults[p] for p ∈ parameters(sys)])

RHs = [10, 25, 40, 55, 65, 70, 75, 80, 85, 90] ./ 100.0
sols = []
for rh in RHs
    p[RH] = rh
    local prob = ODEProblem(sys, u₀, (0, 100.0), p)
    local sol = solve(prob, Rosenbrock23(), abstol=1e-12, reltol=1e-12)
    push!(sols, sol)
end

plot(
    plot(RHs, [sols[i][W][end] * 1e9 for i ∈ 1:length(RHs)], ylim=(0, 50),
        ylabel="H2O (ug/m3)", xlabel="Relative humidity (%)", label=:none),
    plot(RHs, [sols[i][K_aq][end] * 1e6 * mw[K_aq] for i ∈ 1:length(RHs)], #ylim=(0, 1.2),
        ylabel="K+(aq) (ug/m3)", xlabel="Relative humidity (%)", label=:none),
    plot(RHs, [sols[i][NH4_aq][end] * 1e6 * mw[NH4_aq] for i ∈ 1:length(RHs)], #ylim=(0, 3.0),
        ylabel="NH4(aq) (ug/m3)", xlabel="Relative humidity (%)", label=:none),
    plot(RHs, [sols[i][NO3_aq][end] * 1e6 * mw[NO3_aq] for i ∈ 1:length(RHs)], #ylim=(0, 0.25),
        ylabel="NO3(aq) (ug/m3)", xlabel="Relative humidity (%)", label=:none),
)

plot(
    plot(sols[5][t], sols[5][DRH.f_CaNO32],
        ylabel="f_CaNO32", xlabel="time (s)", label=:none),
    plot(sols[5][t], sols[5][CaNO32_s] * 1e6 * mw[CaNO32_s],
        ylabel="CaNO32_s", xlabel="time (s)", label=:none),
    plot(sols[5][t], sols[5][NO3_aq] * 1e6 * mw[NO3_aq],
        ylabel="NO3_aq", xlabel="time (s)", label=:none),
    plot(sols[5][t], sols[5][HNO3_aq] * 1e6 * mw[HNO3_aq],
        ylabel="HNO3_aq", xlabel="time (s)", label=:none),
    plot(sols[5][t], sols[5][HNO3_g] * 1e6 * mw[HNO3_g],
        ylabel="HNO3_g", xlabel="time (s)", label=:none),
    plot(sols[5][t], sols[5][Ca_aq] * 1e6 * mw[Ca_aq],
        ylabel="Ca_aq", xlabel="time (s)", label=:none),
    plot(sols[5][t], sols[5][rxn1.sys.k_rev],
        ylabel="k_rev", xlabel="time (s)", label=:none),
    plot(sols[5][t], sols[5][W] * 1e9,
        ylabel="W (ug/m3)", xlabel="time (s)", label=:none),
)


function plot_mass(RHs, molecs, u₀, title; kwargs...)
    y₀ = zeros(length(RHs))
    y = zeros(length(RHs), length(molecs))
    lab = []
    for (j, (s, x)) in enumerate(molecs)
        y₀ .+= s * u₀[x]
        y[:, j] = s .* [sols[i][x][end] for i ∈ 1:length(RHs)]
        push!(lab, string(x))
    end
    p1 = plot(ylabel="$(title) (mol/m3)", xlabel="Relative humidity (%)"; kwargs...)
    areaplot!(p1, RHs, y, label=permutedims(lab))
    plot!(p1, RHs, y₀, label="u₀ $(title)", color=:black, linewidth=2)
end

plot(
    plot_mass(RHs, K_molecs, u₀, "K"),
    plot_mass(RHs, Ca_molecs, u₀, "Ca"),
    plot_mass(RHs, Mg_molecs, u₀, "Mg"),
    plot_mass(RHs, NH_molecs, u₀, "NH"),
    plot_mass(RHs, Na_molecs, u₀, "Na"),
    plot_mass(RHs, SO4_molecs, u₀, "SO4"),
    plot_mass(RHs, NO3_molecs, u₀, "NO3"),
    plot_mass(RHs, Cl_molecs, u₀, "Cl"),
    plot_mass(RHs, H_molecs, u₀, "H"),
    size=(1000, 800)
)