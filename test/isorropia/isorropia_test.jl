using EarthSciMLBase
using ModelingToolkit, Catalyst, DifferentialEquations, Unitful
using Test
using Plots

include(joinpath(@__DIR__, "../../src/isorropia/isorropia.jl"))

@variables t [unit = u"s", description = "Time"]

model = Isorropia(t)

@test all([ModelingToolkit.check_units(eq) for eq in equations(get_mtk(model))])

sys = structural_simplify(get_mtk(model))

u₀ = ModelingToolkit.get_defaults(sys)
#u₀[RH] = 0.1
#u₀[CaNO32_s] = 2 / 1e6 / mw[CaNO32_s]
prob = ODEProblem(sys, u₀, (0.0, 10000.0), u₀)
@time sol = solve(prob, Rosenbrock23(), abstol=1e-11) # Need low tolerance for mass balance checks to pass.


mw = Dict(Na_aq => 22.989769, SO4_aq => 96.0636, SO4_g => 96.0636, NH3_aq => 17.03052, NH3_g => 17.03052, NO3_aq => 62.0049, Cl_aq => 35.453, 
        Ca_aq => 40.078, K_aq => 39.0983, Mg_aq => 24.305, H_aq => 1.00784, NH4_aq => 18.04,
        K2SO4_s => 174.259, KNO3_s => 101.1032, CaNO32_s => 164.1, HNO3_g => 63.01, HNO3_aq => 63.01) # g/mol


# function plotvars(sol, vars; kwargs...)
#     p1 = plot(; kwargs...)
#     for v ∈ vars
#         plot!(p1, sol[t], sol[v], label=v)
#     end
#     p1
# end

# plot([
# plotvars(sol, [I]; title="I")
# plotvars(sol, [W]; title="W")
# plotvars(sol, states(testsys)[[occursin("aq", string(v)) && (string(v) != "H2O_aq(t)") for v ∈ states(testsys)]];
#     title="Aqueous species")
# plotvars(sol, states(testsys)[[occursin("_s", string(v)) for v ∈ states(testsys)]]; title="Solids")
# plotvars(sol, states(testsys)[[occursin("_g", string(v)) for v ∈ states(testsys)]]; title="Gases")
# plotvars(sol, states(testsys)[[occursin("k_fwd", string(v)) for v ∈ states(testsys)]]; title="Forward rates", yscale=:log10)
# plotvars(sol, states(testsys)[[occursin("k_rev", string(v)) for v ∈ states(testsys)]]; title="Reverse rates", yscale=:log10)
# ]..., size=(1500, 1000))


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

prob = SteadyStateProblem(sys, u₀, u₀)
@time sol = solve(prob, DynamicSS(Rosenbrock23()))
sol[CaNO32_s]  * 1e6 * mw[CaNO32_s]

# Fountoukis and Nenes (2007) Figure 6
u₀ = ModelingToolkit.get_defaults(sys)
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
    prob = SteadyStateProblem(sys, u₀, u₀)
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
