using EarthSciMLBase
using ModelingToolkit, Catalyst, DifferentialEquations, Unitful
using Test
using Plots

include(joinpath(@__DIR__, "../../src/isorropia/isorropia.jl"))
using .ISORROPIA

@variables t [unit = u"s", description = "Time"]

model = Isorropia(t);

sys = structural_simplify(get_mtk(model))

@test all([ModelingToolkit.check_units(eq) for eq in equations(get_mtk(model))])

@unpack Na_aq, SO4_aq, SO4_g, NH3_aq, NH3_g, NO3_aq, Cl_aq, NaCl_s, MgNO32_s, HNO3_g,
    Ca_aq, K_aq, Mg_aq, H_aq, NH4_aq, HCl_g, K2SO4_s, KNO3_s, CaNO32_s, HNO3_g, HNO3_aq,
    KHSO4_s, KCl_s, NH4NO3_s, CaSO4_s, CaCl2_s, MgSO4_s, MgCl2_s, NH4HSO4_s,
    NH42SO4_s, NH43HSO42_s, NH4Cl_s, NaHSO4_s, Na2SO4_s, NaNO3_s, HSO4_aq, HCl_aq,
    RH, metastable, W, I, T, f_CaNO32, rxn1₊k_rev, rxn1₊T₀,
    rxn1₊K⁰, rxn1₊H_group, rxn1₊C_group, T₀₂, c_1, I_one = sys

mw = Dict(Na_aq => 22.989769, SO4_aq => 96.0636, SO4_g => 96.0636, NH3_aq => 17.03052, NH3_g => 17.03052,
    NO3_aq => 62.0049, Cl_aq => 35.453, NaCl_s => 58.44,
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

# eq = []
# for i in eachindex(sol.u)
#     u = [x => sol[x][i] for x in [I, W, states(sys)...]]
#     push!(u, T => sol[T])
#     γ_p = substitute(ModelingToolkit.subs_constants(γ(CaNO32_aqs)), u)
#     a_r = substitute(ModelingToolkit.subs_constants(activity(CaNO32s)), u)
#     a_p = substitute(ModelingToolkit.subs_constants(activity(CaNO32_aqs)), u)
#     keq = sol[rxn1.sys.K_eq][i]
#     @info :a_p => a_p, :γ_p => γ_p, :Ca_aq => sol[Ca_aq][i], :a_r => a_r, :keq => keq, :ratio => a_p / a_r / keq
#     push!(eq, a_p / a_r / keq)
# end

plot(
    plot(sol[t], sol[f_CaNO32],
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
    plot(sol[t], sol[rxn1₊k_rev],
        ylabel="k_rev", xlabel="time (s)", label=:none),
    plot(sol[t], sol[W] * 1e9,
        ylabel="W (ug/m3)", xlabel="time (s)", label=:none),
)

plot([plot(sol[t], sol[ion], ylabel=ion, xlabel="time", label=:none) for ion in ISORROPIA.all_ions]..., size=(1000, 800))
plot([plot(sol[t], sol[gas], ylabel=gas, xlabel="time", label=:none) for gas in ISORROPIA.all_gases]..., size=(1000, 800))
plot([plot(sol[t], sol[solid], ylabel=solid, xlabel="time", label=:none) for solid in ISORROPIA.all_solids]..., size=(1000, 800))

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
        exp(ISORROPIA.logγ₁₂(ISORROPIA.CaNO32_aqs))^(ISORROPIA.CaNO32_aqs.ν_cation + ISORROPIA.CaNO32_aqs.ν_anion)), u)

# First activity calculated using formulas from paper.
aq_activity1 = ca2plus * no3minus^2 * γaq

# Second activity calculated using functions with water conversions.
aq_activity2 = substitute(ModelingToolkit.subs_constants(ISORROPIA.activity(ISORROPIA.CaNO32_aqs)), u)

@test ModelingToolkit.value(aq_activity1) ≈ ModelingToolkit.value(aq_activity2)

keq = rxn1₊K⁰ * exp(-rxn1₊H_group * (rxn1₊T₀ / T - 1) -
                        rxn1₊C_group * (1 + log(rxn1₊T₀ / T) - rxn1₊T₀ / T))
eq_const = substitute(ModelingToolkit.subs_constants(keq), u)

solid_activity = substitute(ModelingToolkit.subs_constants(ISORROPIA.activity(ISORROPIA.CaNO32s)), u)
@test solid_activity ≈ 1.0

# Theoretically the ratio of the activities should be equal to the equilibrium constant
# at the end of the simulation, but it seems like in practice sometimes the concentration on one
# side or the other gets too small and the equilibrium constant blows up.
@test isapprox(ModelingToolkit.value(aq_activity2) /
               ModelingToolkit.value(solid_activity) / ModelingToolkit.value(eq_const), 1.0, atol=1e-2)

# Derivative of activity with respect to concentration should be positive
da_daq = substitute(ModelingToolkit.subs_constants(
        expand_derivatives(Differential(ISORROPIA.Ca_aq)(ISORROPIA.activity(ISORROPIA.CaNO32_aqs)))), [I_one => 1.0, T₀₂ => 273.15, c_1 => 0.005, u...])
@test ModelingToolkit.value(da_daq) > 0.0

# The derivative of our equilibrium ratio (the ratio of our equilibrium expression to one)
# is positive for the aqueous concentration and zero for the solid concentration, because
# the rate of the forward reaction (solid to aqueous) doesn't depend on the concentration of the solid,
# but rate of the reverse reaction (aqueous to solid) does depend on concentration of the aqueous salt.
k_expr = ModelingToolkit.subs_constants(ISORROPIA.activity(ISORROPIA.CaNO32_aqs) / ISORROPIA.activity(ISORROPIA.CaNO32s) / keq)
dk_daq = expand_derivatives(Differential(ISORROPIA.Ca_aq)(k_expr))

@test substitute(dk_daq, u) > 0

dk_ds = expand_derivatives(Differential(ISORROPIA.CaNO32_s)(k_expr))
@test substitute(dk_ds, u₀) == 0

@test ModelingToolkit.get_unit(ISORROPIA.activity(ISORROPIA.CaNO32_aqs) / ISORROPIA.activity(ISORROPIA.CaNO32s) / keq) isa Unitful.FreeUnits{(),NoDims,nothing}


##### Reproducing Figures from Fountoukis and Nenes (2007)

function run_rh_sweep(sys, RHs, ics; mstable=0)
    defaults = ModelingToolkit.get_defaults(sys)
    u₀ = Dict{Any,Float64}([s => 1.0e-15 for s ∈ states(sys)])
    for k ∈ keys(ics)
        u₀[k] = ics[k] / 1e6 / mw[k] # ug/m3 / (1e6 ug/g) / g/mol = mol/m3
    end
    u₀[H_aq] = 2 * u₀[SO4_aq]
    p = Dict{Any,Float64}([p => defaults[p] for p ∈ parameters(sys)])
    p[metastable] = mstable

    sols = []
    for rh in RHs
        p[RH] = rh
        local prob = ODEProblem(sys, u₀, (0, 100.0), p)
        local sol = solve(prob, Rosenbrock23(), abstol=1e-12, reltol=1e-12)
        push!(sols, sol)
    end
    return u₀, sols
end

function plot_rh_sweep(RHs, u₀, sols, plotvars)
    p1 = plot(RHs, [sols[i][W][end] * 1e9 for i ∈ 1:length(RHs)], ylim=(0, 50),
    ylabel="H2O (ug/m3)", xlabel="Relative humidity (%)", label=:none)
    ps = []
    for v in plotvars
        push!(ps, plot(RHs, [sols[i][v][end] * 1e6 * mw[v] for i ∈ 1:length(RHs)],
            ylabel="$v (ug/m3)", xlabel="Relative humidity (%)", label=:none))
    end
    plot(p1, ps...)
end

function plot_mass(RHs, molecs, u₀, sols, title; kwargs...)
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

function plot_all_masses(RHs, u₀, sols)
    plot(
        plot_mass(RHs, K_molecs, u₀, sols, "K"),
        plot_mass(RHs, Ca_molecs, u₀, sols, "Ca"),
        plot_mass(RHs, Mg_molecs, u₀, sols, "Mg"),
        plot_mass(RHs, NH_molecs, u₀, sols, "NH"),
        plot_mass(RHs, Na_molecs, u₀, sols, "Na"),
        plot_mass(RHs, SO4_molecs, u₀, sols, "SO4"),
        plot_mass(RHs, NO3_molecs, u₀, sols, "NO3"),
        plot_mass(RHs, Cl_molecs, u₀, sols, "Cl"),
        plot_mass(RHs, H_molecs, u₀, sols, "H"),
        size=(1000, 800)
    )
end

# Fountoukis and Nenes (2007) Figure 6
RHs = [10, 25, 40, 55, 65, 70, 75, 80, 85, 90] ./ 100.0
ics = Dict([Na_aq => 0, SO4_g => 10, NH3_g => 3.4, HNO3_g => 2, HCl_g => 0,
    Ca_aq => 0.4, K_aq => 0.33, Mg_aq => 1e-20]) # ug/m3
u₀, sols = run_rh_sweep(sys, RHs, ics);
plot_rh_sweep(RHs, u₀, sols, [K_aq, NH4_aq, NO3_aq])
plot_all_masses(RHs, u₀, sols)


# Fountoukis and Nenes (2007) Figure 7
ics = Dict([Na_aq => 3, SO4_g => 3, NH3_g => 0.02, HNO3_g => 2, HCl_g => 3.121,
    Ca_aq => 0.360, K_aq => 0.450, Mg_aq => 0.130]) # ug/m3
u₀, sols = run_rh_sweep(sys, RHs, ics);
plot_rh_sweep(RHs, u₀, sols, [K_aq, NaCl_s, Mg_aq])
plot_all_masses(RHs, u₀, sols)

# Fountoukis and Nenes (2007) Figure 8
ics = Dict([Na_aq => 0.2, SO4_g => 2.0, NH3_g => 8.0, HNO3_g => 12, HCl_g => 0.2,
    Ca_aq => 0.120, K_aq => 0.180, Mg_aq => 0.000]) # ug/m3
u₀, sols = run_rh_sweep(sys, RHs, ics);
plot_rh_sweep(RHs, u₀, sols, [NO3_aq, NH4_aq])
plot_all_masses(RHs, u₀, sols)

# Fountoukis and Nenes (2007) Figure 9
ics = Dict([Na_aq => 0.0, SO4_g => 10.0, NH3_g => 4.250, HNO3_g => 0.145, HCl_g => 0.0,
    Ca_aq => 0.080, K_aq => 0.090, Mg_aq => 0.000]) # ug/m3
u₀1, sols1 = run_rh_sweep(sys, RHs, ics, mstable=0);
#p1 = plot_rh_sweep(RHs, u₀, sols, [K_aq])
u₀2, sols2 = run_rh_sweep(sys, RHs, ics, mstable=1);

p1 = plot(RHs, [sols1[i][W][end] * 1e9 for i ∈ 1:length(RHs)], ylim=(0, 50),
    ylabel="H2O (ug/m3)", xlabel="Relative humidity (%)", label="Stable")
plot!(p1, RHs, [sols2[i][W][end] * 1e9 for i ∈ 1:length(RHs)], ylim=(0, 50), label="Metastable")
ps = []
for v in [K_aq]
    p = plot(RHs, [sols1[i][v][end] * 1e6 * mw[v] for i ∈ 1:length(RHs)],
        ylabel="$v (ug/m3)", xlabel="Relative humidity (%)", label="Stable")
    plot!(p, RHs, [sols2[i][v][end] * 1e6 * mw[v] for i ∈ 1:length(RHs)], label="Metastable")
    push!(ps, p)
end
plot(p1, ps..., size=(600, 300))