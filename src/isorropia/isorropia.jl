module ISORROPIA
using EarthSciMLBase
using ModelingToolkit, Catalyst, Unitful
using IfElse

using Latexify
using NonlinearSolve
using Plots

export Isorropia

include("units.jl")
include("species.jl")
include("aqueous.jl")
include("deliquescence.jl")
include("solid.jl")
include("gas.jl")
include("water.jl")
include("reactions.jl")

# Molar mass in g/mol
mw = Dict(:Na_aq => 22.989769, :SO4_aq => 96.0636, :NH3_aq => 17.03052, :NH3_g => 17.03052,
    :NO3_aq => 62.0049, :Cl_aq => 35.453, :NaCl_s => 58.44,
    :Ca_aq => 40.078, :K_aq => 39.0983, :Mg_aq => 24.305, :H_aq => 1.00784, :NH4_aq => 18.04, :HCl_g => 36.46, :HCl_aq => 36.46,
    :K2SO4_s => :174.259, :KNO3_s => 101.1032, :CaNO32_s => 164.1, :HNO3_g => 63.01, :HNO3_aq => 63.01,
    :KHSO4_s => 136.169, :KCl_s => 74.5513, :NH4NO3_s => 80.043, :NaNO3_s => 84.9947, :NH4Cl_s => 53.491,
    :CaCl2_s => 110.98, :MgCl2_s => 95.211, :NH4HSO4_s => 115.11, :NH42SO4_s => 132.14, :NH43HSO42_s => 247.2485,
    :MgSO4_s => 120.3676, :MgNO32_s => 148.3148, :CaSO4_s => 136.1406, :HSO4_aq => 97.0705, :NaHSO4_s => 120.0603,
    :Na2SO4_s => 142.0421)

# Species containing each molecule for mass balancing
molecs = Dict(
    :K => [(1, :K_aq), (1, :KHSO4_s), (2, :K2SO4_s), (1, :KNO3_s), (1, :KCl_s)],
    :Ca => [(1, :Ca_aq), (1, :CaSO4_s), (1, :CaNO32_s), (1, :CaCl2_s)],
    :Mg => [(1, :Mg_aq), (1, :MgSO4_s), (1, :MgNO32_s), (1, :MgCl2_s)],
    :NH => [(1, :NH4_aq), (1, :NH3_aq), (1, :NH3_g), (1, :NH4HSO4_s),
        (2, :NH42SO4_s), (3, :NH43HSO42_s), (1, :NH4Cl_s), (1, :NH4NO3_s)],
    :Na => [(1, :Na_aq), (1, :NaHSO4_s), (2, :Na2SO4_s), (1, :NaCl_s), (1, :NaNO3_s)],
    :SO4 => [(1, :SO4_aq), (1, :HSO4_aq),
        (1, :KHSO4_s), (1, :NaHSO4_s), (1, :NH4HSO4_s), (1, :CaSO4_s), (1, :Na2SO4_s), (1, :NH42SO4_s),
        (2, :NH43HSO42_s), (1, :K2SO4_s), (1, :MgSO4_s)],
    :NO3 => [(1, :NO3_aq), (1, :HNO3_aq), (1, :HNO3_g), (1, :NH4NO3_s), (1, :NaNO3_s)],
    :Cl => [(1, :Cl_aq), (1, :HCl_aq), (1, :HCl_g), (1, :NH4Cl_s),
        (1, :NaCl_s), (2, :CaCl2_s), (1, :KCl_s), (2, :MgCl2_s)],
    :H => [(1, :H_aq), (1, :HNO3_g), (1, :HCl_g)],
)

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
    sys
    "Molar mass for each species in g/mol"
    mw
    "Species containing each molecule for mass balancing"
    molecs
    "Dictionary of solids"
    solids
    "Dictionary of gases"
    gases
    "Dictionary of aqueous ions"
    ions
    "Dictionary of salts"
    salts

    Isorropia(t; name=:isorropia) = Isorropia(t, :all; name=name) 
    function Isorropia(t, type::Symbol; name=:isorropia)
        sys_types = Dict(
            :all => 1:27, # Use all 27 reactions
            :type1 => [3, 5, 11, 12, 13, 14, 15, 16, 17, 18, 24, 26], # Type 1 reactions (R_1 < 1.0)
        )
        @assert (type ∈ keys(sys_types)) "invalid system type $type, valid options are $(keys(sys_types))"
        Isorropia(t, sys_types[type]; name=name)
    end
    Isorropia(t, rxn_num::Int; name=:isorropia) = Isorropia(t, [rxn_num]; name=name)
    function Isorropia(t, rxn_nums::AbstractVector; name=:isorropia)
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
        active_ions_including_salts = unique(vcat(active_ions, [[s.cation, s.anion] for s ∈ active_salts]...))

        water = Water(t, active_salts)
        ionic_strength = IonicStrength(t, active_ions_including_salts, water.W)
        salt_act = Activities(t, active_salts, (s) -> activity(s, active_salts, salts, ionic_strength.I, water.W))
        ion_act = Activities(t, active_ions, (i) -> activity(i, water.W))
        gas_act = Activities(t, active_gases, activity)
        solid_act = Activities(t, active_solids, activity)
        water_act = Activities(t, [H2Oaq], activity)
        activities = extend(water_act, extend(salt_act, extend(ion_act, extend(gas_act, solid_act))))
        balance = Balance(t, T, active_specs, active_ions_including_salts, gases, solids, ions)
        del, f_del = deliquescence(t, RH, active_salts, salts, ions)
        #out = Output(active_vars)

        # Convert dictionaries of equations into an equation system by combining equations with the same left-hand side.
        eq_dict = Dict()
        syss = []
        for rx ∈ rxns
            sys, ode_eqs_dict = rxn_sys(rx, t, activities, f_del)
            push!(syss, sys)
            for lhs ∈ keys(ode_eqs_dict)
                if lhs in keys(eq_dict)
                    eq_dict[lhs] += ode_eqs_dict[lhs]
                else
                    eq_dict[lhs] = ode_eqs_dict[lhs]
                end
            end
        end
        eqs = [lhs ~ rhs for (lhs, rhs) ∈ pairs(eq_dict)]

        sys = ODESystem(eqs, t, active_vars, [T, RH]; systems=[syss; del], name=name)
        sys = extend(sys, extend(ionic_strength, extend(activities, extend(balance, water))))
        
        filter_present(d, l) = filter(((k,v),) -> !isnothing(findfirst(isequal(v), l)), d)
        new(sys, mw, molecs, 
            filter_present(solids, active_solids), 
            filter_present(gases, active_gases), 
            filter_present(ions, active_ions_including_salts), 
            filter_present(salts, active_salts),
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

    @constants R = 8.31446261815324 [unit = u"m_air^3*Pa/K/mol", description = "Universal gas constant"]
    @constants PaPerAtm = 101325 [unit = u"Pa/atm", description = "Number of pascals per atmosphere"]
    atm2mol(p) = p * PaPerAtm / R / T

    @variables totalK(t) = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of K"]
    @variables totalCa(t) = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of Ca"]
    @variables totalMg(t) = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of Mg"]
    @variables totalNH(t) = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of NH"]
    @variables totalNa(t) = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of Na"]
    @variables totalSO4(t) = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of SO4"]
    @variables totalNO3(t) = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of NO3"]
    @variables totalCl(t) = 1.0e-7 [unit = u"mol/m_air^3", description = "Total concentration of Cl"]
    @variables R_1(t) = 1.0 [description = "Total sulfate ratio from Section 3.1 of Fountoukis and Nenes (2007)"]
    @variables R_2(t) = 1.0 [description = "Crustal species and sodium ratio from Section 3.1 of Fountoukis and Nenes (2007)"]
    @variables R_3(t) = 1.0 [description = "Crustal species ratio from Section 3.1 of Fountoukis and Nenes (2007)"]
    @variables charge(t) = 0 [unit = u"mol/m_air^3", description = "Total ion charge"]
    vs = [totalK, totalCa, totalMg, totalNH, totalNa, totalSO4, totalNO3, totalCl, R_1, R_2, R_3, charge]

    eqs = [
        # Mass balances
        totalK ~ is(i[:K]) + is(sld[:KHSO4]) + 2is(sld[:K2SO4]) + is(sld[:KNO3]) + is(sld[:KCl])
        totalCa ~ is(i[:Ca]) + is(sld[:CaSO4]) + is(sld[:CaNO32]) + is(sld[:CaCl2])
        totalMg ~ is(i[:Mg]) + is(sld[:MgSO4]) + is(sld[:MgNO32]) + is(sld[:MgCl2])
        totalNH ~ is(i[:NH4]) + is(i[:NH3]) + atm2mol(is(g[:NH3])) + is(sld[:NH4HSO4]) +
                  2is(sld[:NH42SO4]) + 3is(sld[:NH43HSO42]) + is(sld[:NH4Cl]) + is(sld[:NH4NO3])
        totalNa ~ is(i[:Na]) + is(sld[:NaHSO4]) + 2is(sld[:Na2SO4]) + is(sld[:NaCl]) + is(sld[:NaNO3])
        totalSO4 ~ is(i[:SO4]) + is(i[:HSO4]) +
                   is(sld[:KHSO4]) + is(sld[:NaHSO4]) + is(sld[:NH4HSO4]) + is(sld[:CaSO4]) + is(sld[:Na2SO4]) + is(sld[:NH42SO4]) +
                   2is(sld[:NH43HSO42]) + is(sld[:K2SO4]) + is(sld[:MgSO4])
        totalNO3 ~ is(i[:NO3]) + is(i[:HNO3]) + atm2mol(is(g[:HNO3])) + is(sld[:NH4NO3]) + is(sld[:NaNO3])
        totalCl ~ is(i[:Cl]) + is(i[:HCl]) + atm2mol(is(g[:HCl])) +
                  is(sld[:NH4Cl]) + is(sld[:NaCl]) + 2is(sld[:CaCl2]) + is(sld[:KCl]) + 2is(sld[:MgCl2])
        # Charge balance
        charge ~ sum([i.m * i.z for i in active_ions_including_salts])

        # Ratios from Section 3.1
        # R_1 ~ (totalNH + totalCa + totalK + totalMg + totalNa) / totalSO4
        # R_2 ~ (totalCa + totalK + totalMg + totalNa) / totalSO4
        # R_3 ~ (totalCa + totalK + totalMg) / totalSO4
    ]
    nonzeros = [e.rhs !== 0 for e ∈ eqs]
    ODESystem(eqs[nonzeros], t, vs, [], name=:balance)
end

"""
Create a system of equations for output variables.
"""
function Output(active_vars)
    names = Symbolics.tosymbol.(active_vars, escape=false)
    ovars = [only(@variables $(n) [unit = ModelingToolkit.get_unit(v), description = "Output for $n"]) for
             (n, v) ∈ zip(names, active_vars)]
    eqs = [ov ~ v for (ov, v) ∈ zip(ovars, active_vars)]
    NonlinearSystem(eqs, ovars, [], name=:out)
end

end