using Aerosol
using ModelingToolkit
using Test

@testitem "initialize Equilibrium Constants" begin
    using ModelingToolkit
    @named eq = ISORROPIA.EquilibriumConstants()
end

@testitem "Compile Isorropia" begin
    using ModelingToolkit
    @named isrpa = Isorropia()
    mtkcompile(isrpa)
end

#@testitem "Run Isorropia" begin
using ModelingToolkit
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using NonlinearSolve, SteadyStateDiffEq

@named isrpa = Isorropia()

sys = mtkcompile(isrpa)

equations(sys)
unknowns(sys)

prob = ODEProblem(sys,
    [
        #sys.TotalNH => 1.0e-8,
        # sys.TotalNa => 1.0e-8,
        # sys.TotalCa => 1.0e-8,
        # sys.TotalK => 1.0e-8,
        # sys.TotalMg => 1.0e-8,
        # sys.TotalCl => 1.0e-8,
        # sys.TotalNO3 => 1.0e-8,
        # sys.TotalSO4 => 5.0e-8,

        sys.aq.NH42SO4.M => 1.0e-8,
    ],
        guesses = [
 sys.aq.NO3.m => 2.58766068993378e-6,
 sys.aq.H.m => 0.003967925913764908,
 sys.aq.Na.m => 0.009506196497264207,
 sys.aq.K.m => 0.00031687425644320324,
 sys.aq.SO4.m => 0.0007586406701139102,
 sys.aq.Ca.m => 2.998719036287425e-7,
 sys.aq.HSO4.m => 0.41395322629634074,
 sys.aq.Mg.m => 3598.450503551135,
 sys.aq.Cl.m => 0.0028059719127294233,
 sys.aq.OH.m => 7196.534148456016,
 sys.aq.MgCl2.M => -0.004072877469889674,
 sys.aq.W => 2.657538735666479e-6,
 sys.aq.NH42SO4.M => 1.0e-8,
 sys.aq.NH4HSO4.loga => 0.3364722366212129,
 sys.aq.MgSO4.loga => 11.608235644774552,
 sys.aq.MgCl2.loga => 50.61605005134875,
 sys.aq.NH42SO4.loga => 0.6418538861723947,
 sys.aq.HNO3.loga => 14.73180128983843,
 sys.aq.Na2SO4.loga => -0.7339691750802004,
 sys.aq.MgNO32.loga => 35.45506712678484,
 sys.aq.K2SO4.loga => -4.135166556742356,
 sys.aq.CaCl2.loga => 27.40787756461434,
 sys.aq.NH43HSO42.loga => 3.4011973816621555
    ],
#    initializealg = ShampineCollocationInit(),
    (0.0, 1.0),
    use_scc = false)

#     collect(zip(unknowns(sys), prob.u0))
# filter(x -> x[2] < 0, collect(zip(unknowns(sys), prob.u0)))

sol = solve(prob, Rosenbrock23())

collect(zip(unknowns(sys), sol.u[1]))
sol[sys.aq.W]
sol[sys.aq.I]

let
    vars = [        sys.aq.NH4.m,sys.aq.Na.m,sys.aq.H.m,        sys.aq.Ca.m,        sys.aq.K.m,
        sys.aq.Mg.m,sys.aq.Cl.m,sys.aq.NO3.m,sys.aq.SO4.m,sys.aq.HSO4.m,sys.aq.OH.m,
        sys.aq.NH3.m,sys.aq.HNO3_aq.m,sys.aq.HCl_aq.m]
        collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
vars = [sys.aq.NH4NO3.loga, sys.aq.NH4Cl.loga, sys.aq.NH4HSO4.loga, sys.aq.NH42SO4.loga, sys.aq.NH43HSO42.loga,
    sys.aq.NaCl.loga, sys.aq.Na2SO4.loga, sys.aq.NaNO3.loga, sys.aq.NaHSO4.loga,
    sys.aq.H2SO4.loga, sys.aq.HCl.loga, sys.aq.HNO3.loga, sys.aq.KHSO4.loga, sys.aq.HHSO4.loga,
    sys.aq.CaNO32.loga, sys.aq.CaCl2.loga, sys.aq.CaSO4.loga,
    sys.aq.KHSO4.loga, sys.aq.K2SO4.loga, sys.aq.KNO3.loga, sys.aq.KCl.loga,
    sys.aq.MgSO4.loga, sys.aq.MgNO32.loga, sys.aq.MgCl2.loga]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
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
vars = [sys.eq.r2.logK_eq, sys.eq.r4.logK_eq, sys.eq.r8.logK_eq, sys.eq.r9.logK_eq,
sys.eq.r10.logK_eq, sys.eq.r14.logK_eq,sys.eq.r15.logK_eq, sys.eq.r19.logK_eq,
sys.eq.r20.logK_eq, sys.eq.r26.logK_eq, sys.eq.r27.logK_eq,
]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [sys.s.NH4, sys.s.Na, sys.s.Ca, sys.s.K, sys.s.Mg, sys.s.Cl, sys.s.NO3,
        sys.s.SO4, sys.s.HSO4]
        collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [sys.g.NH3.p, sys.g.HCl.p, sys.g.HNO3.p, sys.g.H2SO4.p,
        sys.g.NH3.M, sys.g.HCl.M, sys.g.HNO3.M, sys.g.H2SO4.M]
        collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [sys.aq.NH3.m, sys.aq.HCl_aq.m, sys.aq.HNO3_aq.m]
        collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

sol[[sys.aqNH + sys.sNH + sys.g.NH3.M - sys.TotalNH
     sys.aqNa + sys.sNa - sys.TotalNa
     sys.aqCa + sys.sCa - sys.TotalCa
     sys.aqK + sys.sK - sys.TotalK
     sys.aqMg + sys.sMg - sys.TotalMg
     sys.aqCl + sys.sCl + sys.g.HCl.M - sys.TotalCl
     sys.aqNO3 + sys.sNO3 + sys.g.HNO3.M - sys.TotalNO3
     sys.aqSO4 + sys.sSO4 + sys.g.H2SO4.M - sys.TotalSO4]]

iprob = prob.f.initializeprob
sol = solve(iprob)
sol.stats

prob = SteadyStateProblem(iprob.f, iprob.u0, iprob.p)
solve(prob, SSRootfind())

equations(iprob.f.sys)

collect(zip(unknowns(iprob.f.sys), round.(sol.u; sigdigits=2), round.(sol.resid; sigdigits=2)))

using BracketingNonlinearSolve

iprob.p

uspan = (ones(14)*1e-20, ones(14)*10)
iprob2 = IntervalNonlinearProblem{false}(iprob.f, uspan, iprob.p)
solve(iprob2)


     #end
using EarthSciMLBase

allvars = filter(
    x -> occursin("₊M", string(x)), EarthSciMLBase.var2symbol.(unknowns(isrpa)))

# Variables for case 1, sulfate rich free acid.
vars1 = [
    # major
    :s₊NaHSO4₊M,
    :s₊NH4HSO4₊M,
    :s₊KHSO4₊M,
    :s₊CaSO4₊M,

    :aq₊Na₊M,
    :aq₊NH4₊M,
    :aq₊H₊M,

    :aq₊HSO4₊M,
    :aq₊SO4₊M,
    :aq₊NO3₊M,
    :aq₊Cl₊M,
    :aq₊Ca₊M,
    :aq₊K₊M,

    # minor
    :g₊NH3₊M,
    :aq₊NH3₊M,
    :aq₊HNO3_aq₊M,
    :aq₊HCl_aq₊M,
]

@assert length(intersect(allvars, vars1)) == length(vars1)

setdiff(allvars, vars1)

end
##### Reproducing Figures from Fountoukis and Nenes (2007)
# function run_rh_sweep(sys, RHs, ics; mstable = 0)
#     defaults = ModelingToolkit.get_defaults(sys)
#     u₀ = Dict{Any, Float64}([s => 1.0e-15 for s in states(sys)])
#     for k in keys(ics)
#         u₀[k] = ics[k] / 1e6 / mw[k] # ug/m3 / (1e6 ug/g) / g/mol = mol/m3
#     end
#     u₀[H_aq] = 2 * u₀[SO4_aq]
#     p = Dict{Any, Float64}([p => defaults[p] for p in parameters(sys)])
#     p[metastable] = mstable

#     sols = []
#     for rh in RHs
#         p[RH] = rh
#         local prob = ODEProblem(sys, u₀, (0, 100.0), p)
#         local sol = solve(prob, Rosenbrock23(), abstol = 1e-12, reltol = 1e-12)
#         push!(sols, sol)
#     end
#     return u₀, sols
# end

# function plot_rh_sweep(RHs, u₀, sols, plotvars)
#     p1 = plot(
#         RHs,
#         [sols[i][W][end] * 1e9 for i in 1:length(RHs)],
#         ylim = (0, 50),
#         ylabel = "H2O (ug/m3)",
#         xlabel = "Relative humidity (%)",
#         label = :none
#     )
#     ps = []
#     for v in plotvars
#         push!(
#             ps,
#             plot(
#                 RHs,
#                 [sols[i][v][end] * 1e6 * mw[v] for i in 1:length(RHs)],
#                 ylabel = "$v (ug/m3)",
#                 xlabel = "Relative humidity (%)",
#                 label = :none
#             )
#         )
#     end
#     plot(p1, ps...)
# end

# function plot_mass(RHs, molecs, u₀, sols, title; kwargs...)
#     y₀ = zeros(length(RHs))
#     y = zeros(length(RHs), length(molecs))
#     lab = []
#     for (j, (s, x)) in enumerate(molecs)
#         y₀ .+= s * u₀[x]
#         y[:, j] = s .* [sols[i][x][end] for i in 1:length(RHs)]
#         push!(lab, string(x))
#     end
#     p1 = plot(ylabel = "$(title) (mol/m3)", xlabel = "Relative humidity (%)"; kwargs...)
#     areaplot!(p1, RHs, y, label = permutedims(lab))
#     plot!(p1, RHs, y₀, label = "u₀ $(title)", color = :black, linewidth = 2)
# end

# function plot_all_masses(RHs, u₀, sols)
#     plot(
#         plot_mass(RHs, K_molecs, u₀, sols, "K"),
#         plot_mass(RHs, Ca_molecs, u₀, sols, "Ca"),
#         plot_mass(RHs, Mg_molecs, u₀, sols, "Mg"),
#         plot_mass(RHs, NH_molecs, u₀, sols, "NH"),
#         plot_mass(RHs, Na_molecs, u₀, sols, "Na"),
#         plot_mass(RHs, SO4_molecs, u₀, sols, "SO4"),
#         plot_mass(RHs, NO3_molecs, u₀, sols, "NO3"),
#         plot_mass(RHs, Cl_molecs, u₀, sols, "Cl"),
#         plot_mass(RHs, H_molecs, u₀, sols, "H"),
#         size = (1000, 800)
#     )
# end

# # Fountoukis and Nenes (2007) Figure 6
# RHs = [10, 25, 40, 55, 65, 70, 75, 80, 85, 90] ./ 100.0
# ics = Dict([
#     Na_aq => 0,
#     SO4_g => 10,
#     NH3_g => 3.4,
#     HNO3_g => 2,
#     HCl_g => 0,
#     Ca_aq => 0.4,
#     K_aq => 0.33,
#     Mg_aq => 1e-20
# ]) # ug/m3
# u₀, sols = run_rh_sweep(sys, RHs, ics);
# plot_rh_sweep(RHs, u₀, sols, [K_aq, NH4_aq, NO3_aq])
# plot_all_masses(RHs, u₀, sols)

# # Fountoukis and Nenes (2007) Figure 7
# ics = Dict([
#     Na_aq => 3,
#     SO4_g => 3,
#     NH3_g => 0.02,
#     HNO3_g => 2,
#     HCl_g => 3.121,
#     Ca_aq => 0.360,
#     K_aq => 0.450,
#     Mg_aq => 0.130
# ]) # ug/m3
# u₀, sols = run_rh_sweep(sys, RHs, ics);
# plot_rh_sweep(RHs, u₀, sols, [K_aq, NaCl_s, Mg_aq])
# plot_all_masses(RHs, u₀, sols)

# # Fountoukis and Nenes (2007) Figure 8
# ics = Dict([
#     Na_aq => 0.2,
#     SO4_g => 2.0,
#     NH3_g => 8.0,
#     HNO3_g => 12,
#     HCl_g => 0.2,
#     Ca_aq => 0.120,
#     K_aq => 0.180,
#     Mg_aq => 0.000
# ]) # ug/m3
# u₀, sols = run_rh_sweep(sys, RHs, ics);
# plot_rh_sweep(RHs, u₀, sols, [NO3_aq, NH4_aq])
# plot_all_masses(RHs, u₀, sols)

# # Fountoukis and Nenes (2007) Figure 9
# ics = Dict([
#     Na_aq => 0.0,
#     SO4_g => 10.0,
#     NH3_g => 4.250,
#     HNO3_g => 0.145,
#     HCl_g => 0.0,
#     Ca_aq => 0.080,
#     K_aq => 0.090,
#     Mg_aq => 0.000
# ]) # ug/m3
# u₀1, sols1 = run_rh_sweep(sys, RHs, ics, mstable = 0);
# #p1 = plot_rh_sweep(RHs, u₀, sols, [K_aq])
# u₀2, sols2 = run_rh_sweep(sys, RHs, ics, mstable = 1);

# p1 = plot(
#     RHs,
#     [sols1[i][W][end] * 1e9 for i in 1:length(RHs)],
#     ylim = (0, 50),
#     ylabel = "H2O (ug/m3)",
#     xlabel = "Relative humidity (%)",
#     label = "Stable"
# )
# plot!(
#     p1,
#     RHs,
#     [sols2[i][W][end] * 1e9 for i in 1:length(RHs)],
#     ylim = (0, 50),
#     label = "Metastable"
# )
# ps = []
# for v in [K_aq]
#     p = plot(
#         RHs,
#         [sols1[i][v][end] * 1e6 * mw[v] for i in 1:length(RHs)],
#         ylabel = "$v (ug/m3)",
#         xlabel = "Relative humidity (%)",
#         label = "Stable"
#     )
#     plot!(
#         p,
#         RHs,
#         [sols2[i][v][end] * 1e6 * mw[v] for i in 1:length(RHs)],
#         label = "Metastable"
#     )
#     push!(ps, p)
# end
# plot(p1, ps..., size = (600, 300))
