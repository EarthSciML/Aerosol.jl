using Aerosol
using ModelingToolkit
using Test
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using NonlinearSolve

@test string.(ISORROPIA.salt_group(isrpa.aq, :NH4, :M_aq)) == [
     "isrpa₊aq₊NH4NO3₊M_aq(t)"
 "isrpa₊aq₊NH4Cl₊M_aq(t)"
 "isrpa₊aq₊NH4HSO4₊M_aq(t)"
 "isrpa₊aq₊NH42SO4₊M_aq(t)"
 "isrpa₊aq₊NH43HSO42₊M_aq(t)"
]

@test string(sum(ISORROPIA.salt_group(isrpa.aq, :NH4, :M_aq) .* ISORROPIA.salt_group(isrpa.aq, :NH4,
ISORROPIA.salt_group_ν(:NH4)) )) == "isrpa₊aq₊NH42SO4₊ν_cation*isrpa₊aq₊NH42SO4₊M_aq(t) + isrpa₊aq₊NH43HSO42₊ν_cation*isrpa₊aq₊NH43HSO42₊M_aq(t) + isrpa₊aq₊NH4Cl₊ν_cation*isrpa₊aq₊NH4Cl₊M_aq(t) + isrpa₊aq₊NH4HSO4₊ν_cation*isrpa₊aq₊NH4HSO4₊M_aq(t) + isrpa₊aq₊NH4NO3₊ν_cation*isrpa₊aq₊NH4NO3₊M_aq(t)"

@test string(ISORROPIA.resid_moles(isrpa.aq, isrpa.TotalNH, :NH4, :NH43HSO42, [:NH4NO3, :NH4Cl], :M_aq, isrpa.M_zero)) ==
"max(isrpa₊TotalNH(t) / isrpa₊aq₊NH43HSO42₊ν_cation - isrpa₊aq₊NH4Cl₊ν_cation*isrpa₊aq₊NH4Cl₊M_aq(t) - isrpa₊aq₊NH4NO3₊ν_cation*isrpa₊aq₊NH4NO3₊M_aq(t), isrpa₊M_zero)"
@test string(ISORROPIA.resid_moles(isrpa.aq, isrpa.TotalSO4, :SO4, :NH43HSO42, [:CaSO4, :MgSO4], :M_aq, isrpa.M_zero)) ==
"max(isrpa₊TotalSO4(t) / isrpa₊aq₊NH43HSO42₊ν_anion - isrpa₊aq₊CaSO4₊ν_anion*isrpa₊aq₊CaSO4₊M_aq(t) - isrpa₊aq₊MgSO4₊ν_anion*isrpa₊aq₊MgSO4₊M_aq(t), isrpa₊M_zero)"

@named isrpa = Isorropia()

sys = mtkcompile(isrpa)

prob = ODEProblem(sys, [],
    initializealg = BrownFullBasicInit(nlsolve = RobustMultiNewton()),
#    guesses = [
       #sys.aq.NH43HSO42.m_eq => 2.7077033643314414
          # sys.aq.H2SO4.m_eq => 164.92775568117182
          #  sys.aq.HNO3.m_eq => 115.4892998645992
        # sys.aq.HNO3_aq.m_eq => 0.05622990807042508
       #      sys.aq.NH3.m_eq => 5.337486773683763
 #sys.aq.NH3_dissociated.m_eq => 13.482210465994047
       #   sys.aq.HCl_aq.m_eq => 0.0017884729025467545
                 #  sys.Cl_eq => 177.4083678490504
                  # sys.Mg_eq => 134.65917365719395
                 #sys.g.NH3.M => 4.6931978621320305e-8
               # sys.g.HNO3.M => 6.003885040024658e-8
                # sys.g.HCl.M => 1.6040822937392659e-7
             # sys.aq.Na.m_eq => 159.17639267952742
               #    sys.Ca_eq => 12.00641294331137
               #     sys.K_eq => 5.883740758313111
          # sys.aq.K2SO4.m_eq => 0.422536422218256
            # sys.aq.KCl.m_eq => 0.010688599911850177
        #  sys.aq.MgNO32.m_eq => 51.84749703931382
         #  sys.aq.CaCl2.m_eq => 6.293771203936254
        #   sys.aq.KHSO4.m_eq => 4.92221503614245
           # sys.aq.KNO3.m_eq => 0.10576427782229897
       #  sys.aq.NH42SO4.m_eq => 2.0880341178608917
         #   sys.aq.NaCl.m_eq => 0.02229007278269932
       #  sys.aq.NH4HSO4.m_eq => 1.18303213727794
       #   sys.aq.CaNO32.m_eq => 5.648111764614805
           #    sys.aq.H.m_eq => 464.84857632882114
          # sys.aq.NaNO3.m_eq => 0.3923179789196486
       #   sys.aq.NaHSO4.m_eq => 156.11480514341997
       #  sys.aq.Na2SO4.m_eq => 1.3234897422025436
          # sys.aq.MgSO4.m_eq => 0.41864747008975794
          # sys.aq.CaSO4.m_eq => 0.06452997476030989
         #  sys.aq.MgCl2.m_eq => 82.39302914779037
           #   sys.aq.OH.m_eq => 13.482210781433075
                 # sys.SO4_eq => 7.024941091463253
# sys.aq.H2O_dissociated.m_eq => 3.1543902732542984e-7
 #   ],
    (0.0, 10.0), use_scc = false)

iprob = prob.f.initializeprob
solve(iprob)


# f_abs = let
#     logname = [occursin("log", string(var)) for var in unknowns(iprob.f.sys)]
#     function f_abs(u0)
#     u0[.!logname] .= abs.(u0[.!logname])
#     return u0
#     end
# end

# iiprob = NonlinearProblem((u, p) -> iprob.f(abs.(u), p), iprob.u0, iprob.p)
unknowns(iprob.f.sys) .=> iprob.u0
isol = solve(iprob, RobustMultiNewton())
unknowns(iprob.f.sys) .=> round.(isol.resid, sigdigits = 2)
# "sys." .* replace.(replace.(string.(unknowns(iprob.f.sys)), ("₊" => ".",)), ("(t)" => "",)) .=> abs.(isol.u)
guesses = unknowns(iprob.f.sys) .=> abs.(isol.u)

# obs = getproperty.(observed(iprob.f.sys), :lhs)
# obs .=> isol[obs]
# obs_hcl  = getproperty.(filter(s -> occursin("HCl", string(s)), observed(iprob.f.sys)), :lhs)
# obs_hcl .=> isol[obs_hcl]

# prob = ODEProblem(sys, [], #guesses = guesses,
#    initializealg = BrownFullBasicInit(nlsolve = RobustMultiNewton()),
#    (0.0, 10.0), use_scc = false)

sol = solve(prob, Rosenbrock23())

[var => round(val; sigdigits = 2) for (var, val) in zip(unknowns(sys), sol.u[1])]


x = sol[[sys.NH.total, sys.g.NH3.M, sys.aq.NH3.M, sys.aq.NH3_dissociated.M]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

x = sol[[sys.Cl.total, sys.g.HCl.M, sys.aq.HCl.M, sys.aq.HCl_aq.M]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

x = sol[[sys.NO3.total, sys.g.HNO3.M, sys.aq.HNO3.M, sys.aq.HNO3_aq.M]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

x = sol[[sys.SO4.total, sys.aq.HSO4_dissociated.M, sys.aq.H2SO4.M]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0 atol=1e-22

x = sol[[sys.aq.NH3_dissociated.M;
ISORROPIA.salt_group(sys.aq, :NH4, :M) .*
    ISORROPIA.salt_group(sys.aq, :NH4, ISORROPIA.salt_group_ν(:NH4))]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

x = sol[[sys.Na.total;
ISORROPIA.salt_group(sys.aq, :Na, :M) .*
    ISORROPIA.salt_group(sys.aq, :Na, ISORROPIA.salt_group_ν(:Na))]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

x = sol[[sys.Ca.total;
ISORROPIA.salt_group(sys.aq, :Ca, :M) .*
    ISORROPIA.salt_group(sys.aq, :Ca, ISORROPIA.salt_group_ν(:Ca))]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

x = sol[[sys.K.total;
ISORROPIA.salt_group(sys.aq, :K, :M) .*
    ISORROPIA.salt_group(sys.aq, :K, ISORROPIA.salt_group_ν(:K))]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

x = sol[[sys.Mg.total;
ISORROPIA.salt_group(sys.aq, :Mg, :M) .*
    ISORROPIA.salt_group(sys.aq, :Mg, ISORROPIA.salt_group_ν(:Mg))]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

x = sol[[sys.aq.HCl.M;
ISORROPIA.salt_group(sys.aq, :Cl, :M) .*
    ISORROPIA.salt_group(sys.aq, :Cl, ISORROPIA.salt_group_ν(:Cl))]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

ISORROPIA.salt_group(sys.aq, :NO3, :M)
x = sol[[sys.aq.HNO3.M;
ISORROPIA.salt_group(sys.aq, :NO3, :M) .*
    ISORROPIA.salt_group(sys.aq, :NO3, ISORROPIA.salt_group_ν(:NO3))]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

x = sol[[sys.aq.H2SO4.M;
ISORROPIA.salt_group(sys.aq, :SO4, :M) .*
    ISORROPIA.salt_group(sys.aq, :SO4, ISORROPIA.salt_group_ν(:SO4))]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0

x = sol[[sys.aq.HSO4_dissociated.M;
ISORROPIA.salt_group(sys.aq, :HSO4, :M) .*
    ISORROPIA.salt_group(sys.aq, :HSO4, ISORROPIA.salt_group_ν(:HSO4))]][end]
@test all(x .>= 0)
@test x[1] - sum(x[2:end]) ≈ 0.0




let
    vars = [sys.aq.W, sys.aq.I]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

salts = [:CaNO32, :CaCl2, :CaSO4, :KHSO4, :K2SO4, :KNO3, :KCl, :MgSO4, :MgNO32, :MgCl2,
    :NaCl, :Na2SO4, :NaNO3,
    :NH42SO4, :NH4NO3, :NH4Cl, :NH4HSO4, :NaHSO4, :NH43HSO42, :HHSO4, :HNO3, :HCl,
    :H2SO4, :NH3_dissociated,:H2O_dissociated, :HSO4_dissociated]

ions = [:Na, :H, :Ca, :K, :Mg, :OH, :NH3, :HNO3_aq, :HCl_aq]

let
    vars = [reduce(getproperty, [sys, :aq, salt, :loga_eq]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :m_aq]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [sys.aq.F_Ca, sys.aq.F_K, sys.aq.F_Mg, sys.aq.F_Na,
        sys.aq.F_NH4, sys.aq.F_Cl, sys.aq.F_NO3, sys.aq.F_SO4, sys.aq.F_HSO4, sys.aq.F_OH]
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
    vars = [reduce(getproperty, [sys, :aq, salt, :M_aq]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :deliquesced]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :M_salt]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :M_precip]) for salt in salts]
    collect(zip(vars, round.(sol[vars][end]; sigdigits = 2)))
end

let
    vars = [reduce(getproperty, [sys, :aq, salt, :M_dissolved]) for salt in salts]
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
    vars = [reduce(getproperty, [sys, :aq, ion, :m_eq]) for ion in ions]
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


    function f_cat(s)
        x = reduce(getproperty, [sys, :aq, s])
        [x.ν_cation * x.M_aq, x.ν_cation * x.M_precip]
    end
    function f_an(s)
        x = reduce(getproperty, [sys, :aq, s])
        [x.ν_anion * x.M_aq, x.ν_anion * x.M_precip]
    end

    specs = [
        [sys.TotalNH; sys.NH_resid; sys.g.NH3.M; sys.aq.NH3.M;
        vcat(f_cat.([:NH4NO3, :NH4Cl, :NH4HSO4, :NH42SO4, :NH43HSO42])...)],
        [sys.TotalNa;
        vcat(f_cat.([:NaCl, :Na2SO4, :NaNO3, :NaHSO4])...)],
        [sys.TotalCa;
        vcat(f_cat.([:CaNO32, :CaCl2, :CaSO4])...)],
        [sys.TotalK;
        vcat(f_cat.([:KHSO4, :K2SO4, :KNO3, :KCl])...)],
        [sys.TotalMg;
        vcat(f_cat.([:MgSO4, :MgNO32, :MgCl2])...)],
        [sys.TotalCl; sys.g.HCl.M; sys.aq.HCl_aq.M;
        vcat(f_an.([ :NaCl, :KCl, :MgCl2, :CaCl2, :NH4Cl])...)],
        [sys.TotalNO3; sys.g.HNO3.M; sys.aq.HNO3_aq.M;
vcat(f_an.([:NaNO3, :KNO3, :MgNO32, :CaNO32, :NH4NO3])...)],
        [sys.TotalSO4; sys.aq.HHSO4.M_aq;
        vcat(f_an.([:Na2SO4, :K2SO4, :MgSO4, :CaSO4,:NH42SO4, :NH43HSO42,
        :KHSO4, :NaHSO4, :NH4HSO4, :NH43HSO42])...)],
    ]

for i in eachindex(specs)
@info string(specs[i][1])
x = sol[specs[i]][end]
@test all(x .>= 0)
@info string.(specs[i]) .=> x
@info x[1], x[1] - sum(x[2:end])
end

"x"
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
                          guesses = [
        sys.aq.CaNO32.M_salt => 1.906467252413566e-8
       sys.aq.NH43HSO42.m_eq => 2.707703364331441
           sys.aq.H2SO4.m_eq => 164.9277556811717
            sys.aq.HNO3.m_eq => 115.4892998645992
         sys.aq.HNO3_aq.m_eq => 0.056229907643801855
             sys.aq.NH3.m_eq => 5.337486733187579
 sys.aq.NH3_dissociated.m_eq => 13.482210465994047
          sys.aq.HCl_aq.m_eq => 0.0017884728889773864
                   sys.Cl_eq => 177.40836784903678
                 sys.g.NH3.M => 4.693197790930763e-8
                sys.g.HNO3.M => 6.003884948938772e-8
                 sys.g.HCl.M => 1.6040822694034817e-7
              sys.aq.Na.m_eq => 159.17639267952734
                   sys.Ca_eq => 12.006412943311371
                    sys.K_eq => 5.88374075831311
                   sys.Mg_eq => 134.65917365719397
           sys.aq.K2SO4.m_eq => 0.42253642221825605
             sys.aq.KCl.m_eq => 0.010688599911850184
          sys.aq.MgNO32.m_eq => 51.84749703931383
           sys.aq.CaCl2.m_eq => 6.293771203936263
           sys.aq.KHSO4.m_eq => 4.922215036142449
            sys.aq.KNO3.m_eq => 0.10576427782229897
         sys.aq.NH42SO4.m_eq => 2.0880341178608917
            sys.aq.NaCl.m_eq => 0.022290072782699337
         sys.aq.NH4HSO4.m_eq => 1.1830321372779395
          sys.aq.CaNO32.m_eq => 5.648111764614799
               sys.aq.H.m_eq => 464.848576328821
           sys.aq.NaNO3.m_eq => 0.39231797891964865
          sys.aq.NaHSO4.m_eq => 156.1148051434199
          sys.aq.Na2SO4.m_eq => 1.3234897422025436
           sys.aq.MgSO4.m_eq => 0.41864747008975783
           sys.aq.CaSO4.m_eq => 0.06452997476030985
           sys.aq.MgCl2.m_eq => 82.39302914779036
              sys.aq.OH.m_eq => 13.482210781433077
                  sys.SO4_eq => 7.024941091463239
 sys.aq.H2O_dissociated.m_eq => 3.1543902971884954e-7
                      ],
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
