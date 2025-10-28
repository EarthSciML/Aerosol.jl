using Aerosol
using ModelingToolkit
using Test

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
     #        sys.Na_eq => 1e-8,
     sys.TotalNH => 2.3e-7
     sys.TotalNa => 1.0e-8
     sys.TotalCa => 7.5e-10
     sys.TotalK => 3.7e-10
     sys.TotalMg => 8.5e-9
     sys.TotalCl => 1.7e-7
     sys.TotalNO3 => 6.7e-8
     sys.TotalSO4 => 4.4e-10
     # sys.aq.RH => 0.8
     ],
    #          guesses = [
    #            #   sys.aq.NH3.m_aq => 3.1e-7
    #               sys.g.NH3.p => 5.4e-9
    #               sys.g.HNO3.p => 1.4e-9
    #               sys.aq.HNO3_aq.m_aq => 0.0003
    #               sys.g.HCl.p => 3.9e-9
    #               sys.aq.HCl_aq.m_aq => 9.7e-6
    #                  sys.aq.NH3.M => 3.6e-15
    #               sys.aq.HCl_aq.M => 1.1e-13
    #              sys.aq.HNO3_aq.M => 3.5e-12
    #           #  sys.aq.K2SO4.m_eq => 0.42
    #           #    sys.aq.KCl.m_eq => 0.011
    #          #  sys.aq.MgNO32.m_eq => 52.0
    #           #  sys.aq.CaCl2.m_eq => 6.3
    #     #    sys.aq.NH43HSO42.m_eq => 2.7
    #        #     sys.aq.KHSO4.m_eq => 4.9
    #            #  sys.aq.KNO3.m_eq => 0.11
    #           #sys.aq.NH42SO4.m_eq => 2.1
    #          #    sys.aq.NaCl.m_eq => 0.022
    #       #    sys.aq.NH4HSO4.m_eq => 1.2
    #          #  sys.aq.CaNO32.m_eq => 5.6
    #                # sys.aq.H.m_eq => 460.0
    #            # sys.aq.NaNO3.m_eq => 0.39
    #           # sys.aq.NaHSO4.m_eq => 160.0
    #           # sys.aq.Na2SO4.m_eq => 1.3
    #           #  sys.aq.MgSO4.m_eq => 0.42
    #           #  sys.aq.CaSO4.m_eq => 0.065
    #        #     sys.aq.MgCl2.m_eq => 82.0
    #             #   sys.aq.OH.m_eq => 13.0
    # #   sys.aq.H2O_dissociated.W_eq => 6.3e-11
    # #                sys.aq.Na.m_eq => 160.0
    # #                sys.aq.Ca.m_eq => 12.0
    # #                 sys.aq.K.m_eq => 5.9
    # #                sys.aq.Mg.m_eq => 130.0
    #            #   sys.aq.NH3.m_eq => 5.8e-5
    #         #  sys.aq.HNO3_aq.m_eq => 0.056
    #         #   sys.aq.HCl_aq.m_eq => 0.0018
    #           #  sys.aq.H2SO4.m_eq => 160.0
    #            #  sys.aq.HNO3.m_eq => 120.0
    #             #  sys.aq.HCl.m_eq => 180.0
    # #  sys.aq.HSO4_dissociated.m_eq => 7.0
    # #   sys.aq.NH3_dissociated.m_eq => 13.0
    # #   sys.aq.H2O_dissociated.m_eq => 3.2e-7
    #      ],
    initializealg = BrownFullBasicInit(nlsolve = RobustMultiNewton()),
    (0.0, 10.0),
    use_scc = false)

sol = solve(prob, Rosenbrock23())

[var => round(val; sigdigits = 2) for (var, val) in zip(unknowns(sys), sol.u[1])]

let
    vars = [abs(sys.aq.W), sys.aq.W_eq, sys.aq.I]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

salts = [:CaNO32, :CaCl2, :CaSO4, :KHSO4, :K2SO4, :KNO3, :KCl, :MgSO4, :MgNO32, :MgCl2,
    :NaCl, :Na2SO4, :NaNO3,
    :NH42SO4, :NH4NO3, :NH4Cl, :NH4HSO4, :NaHSO4, :NH43HSO42, :HHSO4, :HNO3, :HCl]

ions = [:NH3, :NH3_dissociated, :Na, :Ca, :K, :Mg, :Cl, :NO3, :SO4, :HNO3_aq,
    :HCl_aq, :HSO4,
    :HSO4_dissociated, :H2O_dissociated, :H, :OH]

let
    vars = [reduce(getproperty, [sys, :aq, salt, :loga_eq]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :m_aq]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :m_eq]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :M_eq]) for salt in salts]
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
end

let
    vars = [sys.NH_extra, sys.Na_extra, sys.Ca_extra, sys.K_extra,
        sys.Mg_extra, sys.Cl_extra, sys.NO3_extra, sys.SO4_extra]
    collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [
        sys.TotalNH - sum([sys.aq.NH4NO3.M_aq, sys.aq.NH4Cl.M_aq, sys.aq.NH4HSO4.M_aq,
            2sys.aq.NH42SO4.M_aq, 3sys.aq.NH43HSO42.M_aq]) - sys.g.NH3.M - sys.NH_extra,
        sys.TotalNa - sum([sys.aq.NaCl.M_aq, 2sys.aq.Na2SO4.M_aq, sys.aq.NaNO3.M_aq,
            sys.aq.NaHSO4.M_aq]) - sys.Na_extra,
        sys.TotalCa - sum([sys.aq.CaNO32.M_aq, sys.aq.CaCl2.M_aq, sys.aq.CaSO4.M_aq]) -
        sys.Ca_extra,
        sys.TotalK - sum([sys.aq.KHSO4.M_aq, 2sys.aq.K2SO4.M_aq, sys.aq.KNO3.M_aq,
            sys.aq.KCl.M_aq]) - sys.K_extra,
        sys.TotalMg - sum([sys.aq.MgSO4.M_aq, sys.aq.MgNO32.M_aq, sys.aq.MgCl2.M_aq]) -
        sys.Mg_extra,
        sys.TotalCl - sum([sys.aq.NaCl.M_aq, sys.aq.KCl.M_aq, 2sys.aq.MgCl2.M_aq,
            2sys.aq.CaCl2.M_aq, sys.aq.NH4Cl.M_aq]) - sys.g.HCl.M - sys.Cl_extra,
        sys.TotalNO3 - sum([sys.aq.NaNO3.M_aq, sys.aq.KNO3.M_aq, 2sys.aq.MgNO32.M_aq,
            2sys.aq.CaNO32.M_aq, sys.aq.NH4NO3.M_aq]) - sys.g.HNO3.M - sys.NO3_extra,
        sys.TotalSO4 - sum([sys.aq.Na2SO4.M_aq, sys.aq.K2SO4.M_aq, sys.aq.MgSO4.M_aq,
            sys.aq.CaSO4.M_aq, sys.aq.NH42SO4.M_aq, sys.aq.NH43HSO42.M_aq]) -
        sum([sys.aq.KHSO4.M_eq, sys.aq.NaHSO4.M_eq, sys.aq.NH4HSO4.M_eq,
            sys.aq.NH43HSO42.M_eq]) - sys.SO4_extra
    ]
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
    vars = [sys.aq.NH3_dissociated.M_aq, sys.aq.NH3.M, sys.g.NH3.M,
        sys.aq.NH3_dissociated.m_aq, sys.aq.NH3.m_aq, sys.g.NH3.p]
    collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [sys.aq.HNO3.loga_aq, sys.aq.HNO3.logγ,
        sys.aq.HNO3.M_aq, sys.aq.HNO3_aq.M, sys.g.HNO3.M,
        sys.aq.HNO3.m_aq, sys.aq.HNO3_aq.m_aq, sys.g.HNO3.p]
    collect(zip(vars, round.(sol[vars][1]; sigdigits = 2)))
end

let
    vars = [
        sys.aq.HCl.loga_aq, sys.aq.HCl.logγ, sys.aq.HCl.M_aq, sys.aq.HCl_aq.M, sys.g.HCl.M,
        sys.aq.HCl.m_aq, sys.aq.HCl_aq.m_aq, sys.g.HCl.p]
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

@testset "Cases from paper" begin
    cases = (
        Case = 1:16,
        AerosolType = [
            "Urban (1)", "Urban (2)", "Urban (3)", "Urban (4)",
            "N-u Cont. (1)", "N-u Cont. (2)", "N-u Cont. (3)", "N-u Cont. (4)",
            "Marine (1)", "Marine (2)", "Marine (3)", "Marine (4)",
            "Rem. Cont. (1)", "Rem. Cont. (2)", "Rem. Cont. (3)", "Rem. Cont. (4)"
        ],
        Na = [0.000, 0.023, 0.000, 0.000, 0.200, 0.100, 0.023, 0.023,
            2.000, 1.500, 2.500, 3.000, 0.000, 0.023, 0.100, 0.200],
        H2SO4 = [10.000, 10.000, 15.000, 15.000, 2.000, 4.000, 5.664, 5.664,
            1.000, 1.000, 3.000, 3.000, 10.000, 10.000, 15.000, 15.000],
        NH3 = [3.400, 3.400, 2.000, 2.000, 8.000, 10.000, 12.000, 20.400,
            0.010, 0.010, 0.001, 0.020, 4.250, 3.000, 3.000, 3.000],
        HNO3 = [2.000, 2.000, 10.000, 10.000, 12.000, 7.000, 2.000, 0.611,
            0.300, 1.500, 3.000, 2.000, 0.145, 1.000, 4.000, 8.000],
        HCl = [0.000, 0.037, 0.000, 0.000, 0.200, 0.100, 0.037, 0.037,
            3.121, 2.500, 2.500, 3.121, 0.000, 0.037, 0.100, 0.200],
        Ca2 = [0.400, 0.900, 0.900, 0.400, 0.120, 0.120, 0.120, 0.120,
            0.100, 0.360, 0.500, 0.360, 0.080, 0.080, 0.080, 0.080],
        K = [0.330, 1.000, 1.000, 0.330, 0.180, 0.180, 0.180, 0.180,
            0.100, 0.450, 1.000, 0.450, 0.090, 0.090, 0.090, 0.090],
        Mg2 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.050, 0.050, 0.000,
            0.070, 0.050, 0.050, 0.130, 0.000, 0.000, 0.000, 0.040],
        R1 = [2.14, 2.44, 1.27, 0.89, 23.9, 14.8, 12.4, 20.9,
            9.36, 8.66, 4.86, 5.14, 2.49, 1.78, 1.21, 1.25],
        R2 = [0.18, 0.48, 0.31, 0.12, 0.80, 0.34, 0.18, 0.15,
            9.30, 8.60, 4.86, 5.10, 0.04, 0.05, 0.06, 0.10],
        R3 = [0.18, 0.47, 0.32, 0.12, 0.37, 0.24, 0.17, 0.13,
            0.80, 2.21, 1.31, 0.84, 0.04, 0.04, 0.03, 0.04]
    )

    minval = 0.01

    for i in cases.Case
        @testset "$(cases.AerosolType[i])" begin
            prob = ODEProblem(sys,
                [
                    sys.TotalNa => max(cases.Na[i], minval) / 22.9897693 * 1e-6,
                    sys.TotalSO4 => max(cases.H2SO4[i], minval) / 98.08 * 1e-6,
                    sys.TotalNH => max(cases.NH3[i], minval) / 17.031 * 1e-6,
                    sys.TotalNO3 => max(cases.HNO3[i], minval) / 63.013 * 1e-6,
                    sys.TotalCl => max(cases.HCl[i], minval) / 36.46 * 1e-6,
                    sys.TotalCa => max(cases.Ca2[i], minval) / 40.08 * 1e-6,
                    sys.TotalK => max(cases.K[i], minval) / 39.0983 * 1e-6,
                    sys.TotalMg => max(cases.Mg2[i], minval) / 24.305 * 1e-6
                ],
                #          guesses = [
                #               sys.aq.NH3.m_aq => 3.1e-7
                #               sys.g.NH3.p => 5.4e-9
                #               sys.g.HNO3.p => 1.4e-9
                #               sys.aq.HNO3_aq.m_aq => 0.0003
                #               sys.g.HCl.p => 3.9e-9
                #               sys.aq.HCl_aq.m_aq => 9.7e-6
                #                  sys.aq.NH3.M => 3.6e-15
                #               sys.aq.HCl_aq.M => 1.1e-13
                #              sys.aq.HNO3_aq.M => 3.5e-12
                #             sys.aq.K2SO4.m_eq => 0.42
                #               sys.aq.KCl.m_eq => 0.011
                #            sys.aq.MgNO32.m_eq => 52.0
                #             sys.aq.CaCl2.m_eq => 6.3
                #         sys.aq.NH43HSO42.m_eq => 2.7
                #             sys.aq.KHSO4.m_eq => 4.9
                #              sys.aq.KNO3.m_eq => 0.11
                #           sys.aq.NH42SO4.m_eq => 2.1
                #              sys.aq.NaCl.m_eq => 0.022
                #           sys.aq.NH4HSO4.m_eq => 1.2
                #            sys.aq.CaNO32.m_eq => 5.6
                #                 sys.aq.H.m_eq => 460.0
                #             sys.aq.NaNO3.m_eq => 0.39
                #            sys.aq.NaHSO4.m_eq => 160.0
                #            sys.aq.Na2SO4.m_eq => 1.3
                #             sys.aq.MgSO4.m_eq => 0.42
                #             sys.aq.CaSO4.m_eq => 0.065
                #             sys.aq.MgCl2.m_eq => 82.0
                #                sys.aq.OH.m_eq => 13.0
                #   sys.aq.H2O_dissociated.W_eq => 6.3e-11
                #                sys.aq.Na.m_eq => 160.0
                #                sys.aq.Ca.m_eq => 12.0
                #                 sys.aq.K.m_eq => 5.9
                #                sys.aq.Mg.m_eq => 130.0
                #               sys.aq.NH3.m_eq => 5.8e-5
                #           sys.aq.HNO3_aq.m_eq => 0.056
                #            sys.aq.HCl_aq.m_eq => 0.0018
                #             sys.aq.H2SO4.m_eq => 160.0
                #              sys.aq.HNO3.m_eq => 120.0
                #               sys.aq.HCl.m_eq => 180.0
                #  sys.aq.HSO4_dissociated.m_eq => 7.0
                #   sys.aq.NH3_dissociated.m_eq => 13.0
                #   sys.aq.H2O_dissociated.m_eq => 3.2e-7
                #      ],
                initializealg = BrownFullBasicInit(nlsolve = RobustMultiNewton()),
                (0.0, 10.0),
                use_scc = false)

            sol = solve(prob, Rosenbrock23())

            [var => round(val; sigdigits = 2)
             for (var, val) in zip(unknowns(sys), sol.u[1])]

            iprob = prob.f.initializeprob
            iprob.u0

            ff(u, p) = iprob.f(abs.(u), p)
            iiprob = NonlinearProblem(ff, iprob.u0, iprob.p)
            isol = solve(iiprob, RobustMultiNewton())

            isol.stats

            guesses = unknowns(iprob.f.sys) .=> abs.(isol.u)

            prob = ODEProblem(sys,
                [
                    sys.TotalNa => 0.01 / 22.9897693 * 1e-6,
                    sys.TotalSO4 => 10 / 98.08 * 1e-6,
                    sys.TotalNH => 3.4 / 17.031 * 1e-6,
                    sys.TotalNO3 => 2 / 63.013 * 1e-6,
                    sys.TotalCl => 0.01 / 36.46 * 1e-6,
                    sys.TotalCa => 0.4 / 40.08 * 1e-6,
                    sys.TotalK => 0.33 / 39.0983 * 1e-6,
                    sys.TotalMg => 0.01 / 24.305 * 1e-6
                ],
                guesses = guesses,
                initializealg = BrownFullBasicInit(nlsolve = RobustMultiNewton()),
                (0.0, 10.0),
                use_scc = false)

            sol = solve(prob, Rosenbrock23())

            @test sol.retcode == SciMLBase.ReturnCode.Success
        end
    end
end
