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
        sys.TotalNH => 1e-8
        sys.TotalNa => 1e-8
        sys.TotalCa => 1e-8
        sys.TotalK => 1e-8
        sys.TotalMg => 1e-8
         sys.TotalCl => 1e-8
         sys.TotalNO3 => 1e-8
        sys.TotalSO4 => 8e-8
        sys.aq.RH => 0.3
    ],
         guesses = [
              sys.aq.NH3.a => 6.0e-7
              sys.aq.NH3_dissociated.a => 1.8e-6
              sys.g.HNO3.p => 4.0e-11
              sys.aq.HNO3_aq.a => 0.00085
              sys.aq.HCl_aq.a => 2.6e-5
              sys.aq.H2O_dissociated.a => 5.5e-8
              sys.g.NH3.M => 4.3e-9
              sys.aq.NH4.m => 12.0
              sys.aq.NO3.m => 650.0
                sys.aq.H.m => 0.99
               sys.aq.Na.m => 90.0
                sys.aq.K.m => 3.0
              sys.aq.SO4.m => 5.6
               sys.aq.Ca.m => 16.0
             sys.aq.HSO4.m => 96.0
               sys.aq.Mg.m => 1000.0
               sys.aq.Cl.m => 1400.0
               sys.aq.OH.m => 1.9e-6
  sys.aq.H2O_dissociated.W => 7.2e-12
 sys.aq.HSO4_dissociated.m => 0.99
           sys.aq.CaNO32.m => 7.1
            sys.aq.CaCl2.m => 9.0
            sys.aq.CaSO4.m => 0.03
            sys.aq.K2SO4.m => 0.068
             sys.aq.KNO3.m => 0.015
              sys.aq.KCl.m => 0.00029
            sys.aq.MgSO4.m => 0.0083
           sys.aq.MgNO32.m => 320.0
            sys.aq.MgCl2.m => 690.0
             sys.aq.NaCl.m => 0.0006
           sys.aq.Na2SO4.m => 0.38
            sys.aq.NaNO3.m => 0.056
          sys.aq.NH42SO4.m => 0.75
             sys.aq.HNO3.m => 0.0016
              sys.aq.HCl.m => 1.4e-5
            sys.aq.KHSO4.m => 2.8
          sys.aq.NH4HSO4.m => 0.67
           sys.aq.NaHSO4.m => 89.0
        sys.aq.NH43HSO42.m => 3.4
     ],
    initializealg = BrownFullBasicInit(nlsolve=RobustMultiNewton()),
    (0.0, 10.0),
    use_scc = false)

sol = solve(prob, Rosenbrock23(), abstol = 1e-10, reltol = 1e-10)

[var => round(val; sigdigits=2) for (var, val) in zip(unknowns(sys), sol.u[1])]


let
    vars = [abs(sys.aq.W), sys.aq.I]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

salts = [:CaNO32, :CaCl2, :CaSO4, :KHSO4, :K2SO4, :KNO3, :KCl, :MgSO4, :MgNO32, :MgCl2,
        :NaCl, :Na2SO4, :NaNO3,
        :NH42SO4, :NH4NO3, :NH4Cl, :NH4HSO4, :NaHSO4, :NH43HSO42, :HHSO4, :HNO3, :HCl]

ions = [:NH4, :NH3, :NH3_dissociated, :Na, :Ca, :K, :Mg, :Cl, :NO3, :SO4, :HNO3_aq,
    :HCl_aq, :HSO4,
    :HSO4_dissociated, :H2O_dissociated, :H, :OH]

let
    vars = [reduce(getproperty, [sys, :aq, salt, :loga]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :m]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :M]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end


let
vars = [sys.NH_eq, sys.Na_eq, sys.Ca_eq, sys.K_eq,
    sys.Mg_eq, sys.Cl_eq, sys.NO3_eq, sys.SO4_eq]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
vars = [sys.TotalNH, sys.TotalNa, sys.TotalCa, sys.TotalK,
    sys.TotalMg, sys.TotalCl, sys.TotalNO3, sys.TotalSO4]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
round.(sol[vars][1]; sigdigits = 2)
end

let
vars = [sys.NH_extra, sys.Na_extra, sys.Ca_extra, sys.K_extra,
    sys.Mg_extra, sys.Cl_extra, sys.NO3_extra, sys.SO4_extra]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
vars = [sys.NH_eq + sys.NH_extra - sys.TotalNH,
    sys.Na_eq + sys.Na_extra - sys.TotalNa,
    sys.Ca_eq + sys.Ca_extra - sys.TotalCa,
    sys.K_eq + sys.K_extra - sys.TotalK,
    sys.Mg_eq + sys.Mg_extra - sys.TotalMg,
    sys.Cl_eq + sys.Cl_extra - sys.TotalCl,
    sys.NO3_eq + sys.NO3_extra - sys.TotalNO3,
    sys.SO4_eq + sys.SO4_extra - sys.TotalSO4]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
vars = [reduce(getproperty, [sys, :eq, Symbol(:r, i), :logK_eq]) for i in 1:27]
collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [sys.s.NH4, sys.s.Na, sys.s.Ca, sys.s.K, sys.s.Mg, sys.s.Cl, sys.s.NO3,
        sys.s.SO4, sys.s.HSO4]
        collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, ion, :m]) for ion in ions]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, ion, :M]) for ion in ions]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [sys.g.NH3.p, sys.g.HCl.p, sys.g.HNO3.p,
        sys.g.NH3.M, sys.g.HCl.M, sys.g.HNO3.M]
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
