
#model = Isorropia(t, :all);
#rxn_nums = [10, 11, 12]
rxn_nums = 1:27
#rxn_nums = [10, 11, 12, 13, 21, 26, 27]
#rxn_nums = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
#rxn_nums = 1:27
# rxn 11, 13?
model = Isorropia(t, rxn_nums);

sys = structural_simplify(get_mtk(model))

defaults = ModelingToolkit.get_defaults(sys)
u₀ = Dict{Any,Float64}([s => 1.0e-7 for s ∈ states(sys)])
#u₀[sys.NH3_g] = 0.8 / 1e6 / mw[CaNO32_s]
#u₀[CaNO32_s] = 0.8 / 1e6 / mw[CaNO32_s]
#u₀[Ca_aq] = 0.8 / 1e6 / mw[Ca_aq]
#u₀[NO3_aq] = 10.0 / 1e6 / mw[NO3_aq]

p = Dict{Any,Float64}([p => defaults[p] for p ∈ parameters(sys)])
#p[RH] = 0.30
prob = ODEProblem(sys, u₀, (0.0, 5), p)
# Need low tolerance for mass balance checks to pass.
@time sol = solve(prob, Vern6(), abstol=1e-12, reltol=1e-12;
    callback = PositiveDomain(zeros(length(prob.u0))))

let
    xscale = :none
    yscale = :none
    p1 = plot(title="Solids", xscale=xscale, yscale=yscale, legend=:outertopright)
    for (n, s) in model.solids
        plot!(sol.t[2:end], sol[s.m, 2:end], label=string(n))
    end
    p2 = plot(title="Aqueous Ions", xscale=xscale, yscale=yscale, legend=:outertopright)
    for (n, i) in model.ions
        plot!(sol.t[2:end], sol[i.m, 2:end], label=string(n))
    end
    p3 = plot(title="Gases", xscale=xscale, yscale=yscale)
    for (n, g) in model.gases
        plot!(sol.t[2:end
        ], sol[g.p, 2:end], label=string(n))
    end
    p4 = plot(title="Reaction Rates", xscale=xscale, legend=:outertopright)
    for i ∈ rxn_nums
        r = Symbol(:rxn, i)
        y = eval(:(sol[sys.$r.rate]))
        plot!(sol.t[2:end], y[2:end], label="$r.rate")
        #y2 = eval(:(sol[sys.$r.rawrate]))
        #plot!(sol.t[2:end], y2[2:end], label="$r.rawrate")
    end
    plot(p1, p2, p3, p4, size=(1000, 800))
end

# plot(
#     plot(sol.t, sol[sys.a_K2SO4_aqs], label=sys.a_K2SO4_aqs, yscale=:log10),
#     plot(sol.t, sol[sys.a_K2SO4_s], label=sys.a_K2SO4_s, yscale=:log10),
# #    begin
# #        plot(sol.t, sol[sys.a_NH3_aq] ./ sol[sys.a_NH3_g], label="a ratio")
#         #plot!(sol.t, sol[sys.a_CaNO32_s] ./ sol[sys.a_CaNO32_aqs], yscale=:log10)
#         plot(sol.t, sol[sys.rxn4.K_eq], label="K_eq", yscale=:log10),
# #    end,
# )

xx = 100000
let
    ps = []
    for i ∈ rxn_nums
        r = Symbol(:rxn, i)
        y = eval(:(sol[sys.$r.rate, end-xx:end]))
        p = plot(sol.t[end-xx:end], y, label="$r.rate")
        #y = eval(:(sol[sys.$r.rawrate]))
        #@info y[1:10]
        #p = plot!(sol.t, y, alpha=1, label="$r.rawrate")
        push!(ps, p)
    end
    plot(ps..., size=(1400, 600))
end

sol[sys.rxn16.rawrate]
plot(sol.t, sol[sys.rxn10.rawrate])

plot(
    begin
        plot(sol.t, sol[sys.a_K2SO4_aqs] ./ sol[sys.a_K2SO4_s], label="a ratio", ylim=(0.01, 0.02))
        #plot!(sol.t, sol[sys.a_CaNO32_s] ./ sol[sys.a_CaNO32_aqs], yscale=:log10)
        plot!(sol.t, sol[sys.rxn4.K_eq], label="K_eq")
    end,
    plot(sol.t, sol[sys.a_K2SO4_aqs] ./ sol[sys.a_K2SO4_s] - sol[sys.rxn4.K_eq], 
        label="a ratio - k_eq", ylim=(-1e-6, 1e-6)),
    scatter(sol.t, sol[sys.rxn4.rate], label="rate", ylim=(-2e-6, 2e-6)),
    plot(sol.t, sol[sys.rxn4.rawrate], label="rawrate"),
)





yy = sol[sys.a_NH3_aq] ./ sol[sys.a_NH3_g] - sol[sys.rxn12.K_eq]

softplus(x) = log(1 + exp(x-10))
trans(x) = sign(x) * softplus(abs(x))
plot(
scatter(sol.t, sol[sys.rxn4.rawrate], label="rawrate"),
plot(sol.t, trans.(sol[sys.rxn4.rawrate]), label="rawrate"),
)
xx = 

@test all([ModelingToolkit.check_units(eq) for eq in equations(get_mtk(model))])


let
    ps = []
    for i ∈ rxn_nums
        r = Symbol(:rxn, i)
        y = eval(:(sol[sys.$r.rate]))
        p = plot(sol.t[2:end], y[2:end], xscale=:log10, label="$r.rate")
        push!(ps, p)
    end
    plot(ps..., size=(1000, 800))
end

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