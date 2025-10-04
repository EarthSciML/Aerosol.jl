module ISORROPIA
using EarthSciMLBase
using ModelingToolkit, Catalyst, DynamicQuantities
using IfElse

using Latexify
using NonlinearSolve
using Plots

export Isorropia

include("aqueous.jl")
include("deliquescence.jl")
include("solid.jl")
include("gas.jl")
include("equilibria.jl")

@mtkmodel Isorropia begin
    @components begin
        aq = Aqueous()
        s = Solids()
        g = Gases()
        eq = EquilibriumConstants()
    end
    @parameters begin
        T, [unit = u"K", description = "Temperature"]
        RH, [unit = u"1", description = "Relative humidity (0 to 1)"]
    end
    @variables begin
        TotalNH(t), [unit = u"mol/m^3", description = "Total concentration of NH3 + NH4+"]
        TotalNa(t), [unit = u"mol/m^3", description = "Total concentration of Na+"]
        TotalCa(t), [unit = u"mol/m^3", description = "Total concentration of Ca2+"]
        TotalK(t), [unit = u"mol/m^3", description = "Total concentration of K+"]
        TotalMg(t), [unit = u"mol/m^3", description = "Total concentration of Mg2+"]
        TotalCl(t), [unit = u"mol/m^3", description = "Total concentration of Cl-"]
        TotalNO3(t), [unit = u"mol/m^3", description = "Total concentration of NO3-"]
        TotalSO4(t), [unit = u"mol/m^3", description = "Total concentration of SO4 2-"]
    end
    @equations begin
        T ~ aq.T
        T ~ eq.T
        T ~ g.T
        aq.W ~ s.W
        RH ~ aq.RH

        # Reactions based on information in Table 2 of Fountoukis and Nenes (2007).
        # The left-hand side of the reaction equation is the equilibrium constant and the
        # right-hand side is the ratio of the product and reactant activities.
        eq.r1.K_eq * eq.k1_unit ~ aq.a_CaNO32
        eq.r2.K_eq * eq.k2_unit ~ aq.a_CaCl2
        # eq.r3.K_eq * eq.k3_unit ~ aq.a_CaSO4 * RH^2
        aq.CaSO4.M ~ 0.0 # From Table 4 footnote a, CaSO4 has an activity coefficient of zero and therefore 0 concentration.
        eq.r4.K_eq * eq.k4_unit ~ aq.a_K2SO4
        eq.r5.K_eq * eq.k5_unit ~ aq.a_KHSO4
        eq.r6.K_eq * eq.k6_unit ~ aq.a_KNO3
        eq.r7.K_eq * eq.k7_unit ~ aq.a_KCl
        eq.r8.K_eq * eq.k8_unit ~ aq.a_MgSO4
        eq.r9.K_eq * eq.k9_unit ~ aq.a_MgNO32
        eq.r10.K_eq * eq.k10_unit ~ aq.a_MgCl2
        eq.r11.K_eq * eq.k11_unit ~ aq.H.a * aq.SO4.a / aq.HSO4.a
        eq.r12.K_eq * eq.k12_unit ~ aq.NH3.a / g.NH3.p
        eq.r13.K_eq * eq.k13_unit ~ aq.NH4.a * aq.OH.a / aq.NH3.a / RH
        eq.r14.K_eq * eq.k14_unit ~ aq.H.a * aq.NO3.a / g.HNO3.p # K1
        eq.r15.K_eq * eq.k15_unit ~ aq.a_HNO3 / g.HNO3.p # K1a
        eq.r14.K_eq / eq.r15.K_eq ~ aq.H.a * aq.NO3.a / aq.a_HNO3 # K1b, from Table 2 footnote ♠
        eq.r16.K_eq * eq.k16_unit ~ aq.a_HCl / g.HCl.p # K2
        eq.r17.K_eq * eq.k17_unit ~ aq.HCl_aq.a / g.HCl.p # K2a
        eq.r16.K_eq / eq.r17.K_eq * eq.k2b_unit ~ aq.H.a * aq.Cl.a / aq.HCl_aq.a # K2b, from Table 2 footnote ♦
        eq.r18.K_eq * eq.k18_unit ~ aq.H.a * aq.OH.a / RH
        eq.r19.K_eq * eq.k19_unit ~ aq.a_Na2SO4
        eq.r20.K_eq * eq.k20_unit ~ aq.a_NH42SO4
        eq.r21.K_eq * eq.k21_unit ~ g.NH3.p * g.HCl.p
        eq.r22.K_eq * eq.k22_unit ~ aq.a_NaNO3
        eq.r23.K_eq * eq.k23_unit ~ aq.a_NaCl
        eq.r24.K_eq * eq.k24_unit ~ aq.a_NaHSO4
        eq.r25.K_eq * eq.k25_unit ~ g.NH3.p * g.HNO3.p
        eq.r26.K_eq * eq.k26_unit ~ aq.a_NH4HSO4
        eq.r27.K_eq * eq.k27_unit ~ aq.a_NH43HSO42

        # Mass Balance
        TotalNH ~ aq.NH4.M + aq.NH3.M + g.NH3.M + s.NH4Cl.M + s.NH4NO3.M + 2s.NH42SO4.M +
            s.NH4HSO4.M + 3s.NH43HSO42.M
        TotalNa ~ aq.Na.M + s.NaCl.M + s.NaNO3.M + 2s.Na2SO4.M + s.NaHSO4.M
        TotalCa ~ aq.Ca.M + s.CaCl2.M + s.CaNO32.M + s.CaSO4.M
        TotalK ~ aq.K.M + s.KCl.M + s.KNO3.M + s.KHSO4.M + 2s.K2SO4.M
        TotalMg ~ aq.Mg.M + s.MgCl2.M + s.MgNO32.M + s.MgSO4.M
        TotalCl ~ aq.Cl.M + aq.HCl.M + g.HCl.M + 2s.CaCl2.M + s.KCl.M + s.NH4Cl.M +
            s.NaCl.M + 2s.MgCl2.M
        TotalNO3 ~ aq.NO3.M + g.HNO3.M + 2s.CaNO32.M + s.KNO3.M + s.NH4NO3.M + s.NaNO3.M +
            2s.MgNO32.M
        TotalSO4 ~ aq.SO4.M + g.H2SO4.M + 2s.NH42SO4.M + s.NH4HSO4.M + 2s.NH43HSO42.M +
            s.NaHSO4.M + s.Na2SO4.M + s.KHSO4.M + 2s.K2SO4.M + s.CaSO4.M + s.MgSO4.M

        D(TotalNH) ~ 0.0
        D(TotalNa) ~ 0.0
        D(TotalCa) ~ 0.0
        D(TotalK) ~ 0.0
        D(TotalMg) ~ 0.0
        D(TotalCl) ~ 0.0
        D(TotalNO3) ~ 0.0
        D(TotalSO4) ~ 0.0
    end
end

@named isrpa = Isorropia()

structural_simplify(isrpa)
sys = structural_simplify(isrpa, fully_determined=false)

unknowns(sys)
equations(sys)

isrpa.T

equations(isrpa)


# Molar mass in g/mol
mw = Dict(
    :Na_aq => 22.989769,
    :SO4_aq => 96.0636,
    :NH3_aq => 17.03052,
    :NH3_g => 17.03052,
    :NO3_aq => 62.0049,
    :Cl_aq => 35.453,
    :NaCl_s => 58.44,
    :Ca_aq => 40.078,
    :K_aq => 39.0983,
    :Mg_aq => 24.305,
    :H_aq => 1.00784,
    :NH4_aq => 18.04,
    :HCl_g => 36.46,
    :HCl_aq => 36.46,
    :K2SO4_s => :174.259,
    :KNO3_s => 101.1032,
    :CaNO32_s => 164.1,
    :HNO3_g => 63.01,
    :HNO3_aq => 63.01,
    :KHSO4_s => 136.169,
    :KCl_s => 74.5513,
    :NH4NO3_s => 80.043,
    :NaNO3_s => 84.9947,
    :NH4Cl_s => 53.491,
    :CaCl2_s => 110.98,
    :MgCl2_s => 95.211,
    :NH4HSO4_s => 115.11,
    :NH42SO4_s => 132.14,
    :NH43HSO42_s => 247.2485,
    :MgSO4_s => 120.3676,
    :MgNO32_s => 148.3148,
    :CaSO4_s => 136.1406,
    :HSO4_aq => 97.0705,
    :NaHSO4_s => 120.0603,
    :Na2SO4_s => 142.0421
)

# Species containing each molecule for mass balancing
molecs = Dict(
    :K => [(1, :K_aq), (1, :KHSO4_s), (2, :K2SO4_s), (1, :KNO3_s), (1, :KCl_s)],
    :Ca => [(1, :Ca_aq), (1, :CaSO4_s), (1, :CaNO32_s), (1, :CaCl2_s)],
    :Mg => [(1, :Mg_aq), (1, :MgSO4_s), (1, :MgNO32_s), (1, :MgCl2_s)],
    :NH => [
        (1, :NH4_aq),
        (1, :NH3_aq),
        (1, :NH3_g),
        (1, :NH4HSO4_s),
        (2, :NH42SO4_s),
        (3, :NH43HSO42_s),
        (1, :NH4Cl_s),
        (1, :NH4NO3_s)
    ],
    :Na => [(1, :Na_aq), (1, :NaHSO4_s), (2, :Na2SO4_s), (1, :NaCl_s), (1, :NaNO3_s)],
    :SO4 => [
        (1, :SO4_aq),
        (1, :HSO4_aq),
        (1, :KHSO4_s),
        (1, :NaHSO4_s),
        (1, :NH4HSO4_s),
        (1, :CaSO4_s),
        (1, :Na2SO4_s),
        (1, :NH42SO4_s),
        (2, :NH43HSO42_s),
        (1, :K2SO4_s),
        (1, :MgSO4_s)
    ],
    :NO3 => [(1, :NO3_aq), (1, :HNO3_aq), (1, :HNO3_g), (1, :NH4NO3_s), (1, :NaNO3_s)],
    :Cl => [
        (1, :Cl_aq),
        (1, :HCl_aq),
        (1, :HCl_g),
        (1, :NH4Cl_s),
        (1, :NaCl_s),
        (2, :CaCl2_s),
        (1, :KCl_s),
        (2, :MgCl2_s)
    ],
    :H => [(1, :H_aq), (1, :HNO3_g), (1, :HCl_g)]
)

abstract type EarthSciMLODESystem end

"""
    Isorropia(t)
    Isorropia(t, 1:3) # Only include the first three of the 27 reactions.
    Isorropia(t, 1) # Only include the first reaction.
    Isorropia(t, :all) # Choose a pre-specified set of reactions. Current options are `:all` and `:type1` (referring to first aerosol type in Fountoukis and Nenes (2007) Table 3).

An implementation of ISORROPIA II, a model for the thermodynamic equilibrium of gas-aerosol interactions, as described in:

> Fountoukis, C. and Nenes, A., 2007. ISORROPIA II: a computationally efficient thermodynamic equilibrium model for K+–Ca 2+–Mg 2+–NH 4+–Na+–SO 4 2−–NO 3−–Cl−–H 2 O aerosols. Atmospheric Chemistry and Physics, 7(17), pp.4639-4659.
"""
struct Isorropia <: EarthSciMLODESystem
    "ModelingToolkit ODESystem"
    sys::Any
    "Molar mass for each species in g/mol"
    mw::Any
    "Species containing each molecule for mass balancing"
    molecs::Any
    "Dictionary of solids"
    solids::Any
    "Dictionary of gases"
    gases::Any
    "Dictionary of aqueous ions"
    ions::Any
    "Dictionary of salts"
    salts::Any

    Isorropia(t; name = :isorropia) = Isorropia(t, :all; name = name)
    function Isorropia(t, type::Symbol; name = :isorropia)
        sys_types = Dict(
            :all => 1:27, # Use all 27 reactions
            :type1 => [3, 5, 11, 12, 13, 14, 15, 16, 17, 18, 24, 26] # Type 1 reactions (R_1 < 1.0)
        )
        @assert (type ∈ keys(sys_types)) "invalid system type $type, valid options are $(keys(sys_types))"
        Isorropia(t, sys_types[type]; name = name)
    end
    Isorropia(t, rxn_num::Int; name = :isorropia) = Isorropia(t, [rxn_num]; name = name)
    function Isorropia(t, rxn_nums::AbstractVector; name = :isorropia)
        ions = generate_ions(t)
        salts = generate_salts(ions)
        gases = generate_gases(t)
        solids = generate_solids(t)
        H2Oaq = H2O(RH)

        rxns = generate_reactions(salts, gases, solids, ions, H2Oaq)[rxn_nums]
        active_specs = active_species(rxns)
        active_vars = unique(vcat(vars.(active_specs)...))
        active_ions = intersect(active_specs, values(ions))
        active_salts = intersect(active_specs, values(salts))
        active_solids = intersect(active_specs, values(solids))
        active_gases = intersect(active_specs, values(gases))
        active_ions_including_salts = unique(vcat(
            active_ions, [[s.cation, s.anion] for s in active_salts]...))

        water = Water(t, active_salts)
        ionic_strength = IonicStrength(t, active_ions_including_salts, water.W)
        salt_act = Activities(
            t,
            active_salts,
            (s) -> activity(s, active_salts, salts, ionic_strength.I, water.W)
        )
        ion_act = Activities(t, active_ions, (i) -> activity(i, water.W))
        gas_act = Activities(t, active_gases, activity)
        solid_act = Activities(t, active_solids, activity)
        water_act = Activities(t, [H2Oaq], activity)
        activities = extend(
            water_act, extend(salt_act, extend(ion_act, extend(gas_act, solid_act))))
        balance = Balance(
            t, T, active_specs, active_ions_including_salts, gases, solids, ions)
        del, f_del = deliquescence(t, RH, active_salts, salts, ions)
        #out = Output(active_vars)

        # Convert dictionaries of equations into an equation system by combining equations with the same left-hand side.
        eq_dict = Dict()
        syss = []
        for rx in rxns
            sys, ode_eqs_dict = rxn_sys(rx, t, activities, f_del)
            push!(syss, sys)
            for lhs in keys(ode_eqs_dict)
                if lhs in keys(eq_dict)
                    eq_dict[lhs] += ode_eqs_dict[lhs]
                else
                    eq_dict[lhs] = ode_eqs_dict[lhs]
                end
            end
        end
        eqs = [lhs ~ rhs for (lhs, rhs) in pairs(eq_dict)]

        sys = ODESystem(eqs, t, active_vars, [T, RH]; systems = [syss; del], name = name)
        sys = extend(
            sys, extend(ionic_strength, extend(activities, extend(balance, water))))

        filter_present(d, l) = filter(((k, v),) -> !isnothing(findfirst(isequal(v), l)), d)
        new(
            sys,
            mw,
            molecs,
            filter_present(solids, active_solids),
            filter_present(gases, active_gases),
            filter_present(ions, active_ions_including_salts),
            filter_present(salts, active_salts)
        )
    end
end

"""
Create a system of equations for mass and charge balances,
where slt, g, sld, and i are dictionaries of gases, solids, and ions, respectively.
"""
function Balance(t, T, active_specs, active_ions_including_salts, g, sld, i)
    # Return the given species, or zero if the species is not in the list.
    is(v) = ((v ∈ active_specs) || (v ∈ active_ions_including_salts)) ? only(vars(v)) : 0

    @constants R = 8.31446261815324 [
        unit = u"m_air^3*Pa/K/mol",
        description = "Universal gas constant"
    ]
    @constants PaPerAtm = 101325 [
        unit = u"Pa/Constants.atm",
        description = "Number of pascals per atmosphere"
    ]
    atm2mol(p) = p * PaPerAtm / R / T

    @variables totalK(t) = 1.0e-7 [
        unit = u"mol/m_air^3",
        description = "Total concentration of K"
    ]
    @variables totalCa(t) = 1.0e-7 [
        unit = u"mol/m_air^3",
        description = "Total concentration of Ca"
    ]
    @variables totalMg(t) = 1.0e-7 [
        unit = u"mol/m_air^3",
        description = "Total concentration of Mg"
    ]
    @variables totalNH(t) = 1.0e-7 [
        unit = u"mol/m_air^3",
        description = "Total concentration of NH"
    ]
    @variables totalNa(t) = 1.0e-7 [
        unit = u"mol/m_air^3",
        description = "Total concentration of Na"
    ]
    @variables totalSO4(t) = 1.0e-7 [
        unit = u"mol/m_air^3",
        description = "Total concentration of SO4"
    ]
    @variables totalNO3(t) = 1.0e-7 [
        unit = u"mol/m_air^3",
        description = "Total concentration of NO3"
    ]
    @variables totalCl(t) = 1.0e-7 [
        unit = u"mol/m_air^3",
        description = "Total concentration of Cl"
    ]
    @variables R_1(t) = 1.0 [
        description = "Total sulfate ratio from Section 3.1 of Fountoukis and Nenes (2007)",
    ]
    @variables R_2(t) = 1.0 [
        description = "Crustal species and sodium ratio from Section 3.1 of Fountoukis and Nenes (2007)",
    ]
    @variables R_3(t) = 1.0 [
        description = "Crustal species ratio from Section 3.1 of Fountoukis and Nenes (2007)",
    ]
    @variables charge(t) = 0 [unit = u"mol/m_air^3", description = "Total ion charge"]
    vs = [
        totalK,
        totalCa,
        totalMg,
        totalNH,
        totalNa,
        totalSO4,
        totalNO3,
        totalCl,
        R_1,
        R_2,
        R_3,
        charge
    ]

    eqs = [
           # Mass balances
           totalK ~
           is(i[:K]) + is(sld[:KHSO4]) + 2is(sld[:K2SO4]) + is(sld[:KNO3]) + is(sld[:KCl])
           totalCa ~ is(i[:Ca]) + is(sld[:CaSO4]) + is(sld[:CaNO32]) + is(sld[:CaCl2])
           totalMg ~ is(i[:Mg]) + is(sld[:MgSO4]) + is(sld[:MgNO32]) + is(sld[:MgCl2])
           totalNH ~
           is(i[:NH4]) +
           is(i[:NH3]) +
           atm2mol(is(g[:NH3])) +
           is(sld[:NH4HSO4]) +
           2is(sld[:NH42SO4]) +
           3is(sld[:NH43HSO42]) +
           is(sld[:NH4Cl]) +
           is(sld[:NH4NO3])
           totalNa ~
           is(i[:Na]) +
           is(sld[:NaHSO4]) +
           2is(sld[:Na2SO4]) +
           is(sld[:NaCl]) +
           is(sld[:NaNO3])
           totalSO4 ~
           is(i[:SO4]) +
           is(i[:HSO4]) +
           is(sld[:KHSO4]) +
           is(sld[:NaHSO4]) +
           is(sld[:NH4HSO4]) +
           is(sld[:CaSO4]) +
           is(sld[:Na2SO4]) +
           is(sld[:NH42SO4]) +
           2is(sld[:NH43HSO42]) +
           is(sld[:K2SO4]) +
           is(sld[:MgSO4])
           totalNO3 ~
           is(i[:NO3]) +
           is(i[:HNO3]) +
           atm2mol(is(g[:HNO3])) +
           is(sld[:NH4NO3]) +
           is(sld[:NaNO3])
           totalCl ~
           is(i[:Cl]) +
           is(i[:HCl]) +
           atm2mol(is(g[:HCl])) +
           is(sld[:NH4Cl]) +
           is(sld[:NaCl]) +
           2is(sld[:CaCl2]) +
           is(sld[:KCl]) +
           2is(sld[:MgCl2])
           # Charge balance
           charge ~ sum([i.m * i.z for i in active_ions_including_salts])

           # Ratios from Section 3.1
           # R_1 ~ (totalNH + totalCa + totalK + totalMg + totalNa) / totalSO4
           # R_2 ~ (totalCa + totalK + totalMg + totalNa) / totalSO4
           # R_3 ~ (totalCa + totalK + totalMg) / totalSO4
           ]
    nonzeros = [e.rhs !== 0 for e in eqs]
    ODESystem(eqs[nonzeros], t, vs, [], name = :balance)
end

"""
Create a system of equations for output variables.
"""
function Output(active_vars)
    names = Symbolics.tosymbol.(active_vars, escape = false)
    ovars = [only(
                 @variables $(n) [
                 unit = ModelingToolkit.get_unit(v),
                 description = "Output for $n"
             ]
             ) for (n, v) in zip(names, active_vars)]
    eqs = [ov ~ v for (ov, v) in zip(ovars, active_vars)]
    NonlinearSystem(eqs, ovars, [], name = :out)
end

end


spc = [
:CaNO32s
:CaNO32aq
:CaCl2s
:CaCl2aq
:CaSO4s
:CaSO4aq
:K2SO4s
:K2SO4aq
:KHSO4s
:KHSO4aq
:KNO3s
:KNO3aq
:KCls
:KClaq
:MgSO4s
:MgSO4aq
:MgNO32s
:MgNO32aq
:MgCl2s
:MgCl2aq
:HSO4s
:Haq
:SO4aq
:NH3g
:NH3aq
:NH3aq
:NH4aq
:OHaq
:HNO3g
:HNO3aq
:HClg
:Haq
:Claq
:HClg
:HClaq
:Haq
:OHaq
:Na2SO4s
:Na2SO4aq
:NH4Cls
:NH3g
:HClg
:NaNO3s
:NaNO3aq
:NaCls
:NaClaq
:NaHSO4s
:NaHSO4aq
:NH4NO3s
:NH3g
:HNO3g
:NH4HSO4s
:NH4HSO4aq
:NH43HSO42s
:NH43HSO42aq
]

unique(spc)
