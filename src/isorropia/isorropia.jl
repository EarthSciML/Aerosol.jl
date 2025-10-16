module ISORROPIA

using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities

export Isorropia

include("aqueous.jl")
#include("deliquescence.jl")
include("solid.jl")
include("gas.jl")
include("equilibria.jl")

@mtkmodel Isorropia begin
    @parameters begin
        T = 293.15, [unit = u"K", description = "Temperature"]
        RH = 0.3, [unit = u"1", description = "Relative humidity (0 to 1)"]

        m_one = 1.0, [unit = u"mol/kg"]
        M_one = 1.0, [unit = u"mol/m^3"]
        p_one = 1.0, [unit = u"Constants.atm"]
    end
    @constants begin
        M_zero = 0.0, [unit = u"mol/m^3"]
    end
    @components begin
        aq = Aqueous(T = T, RH = RH)
        s = Solids()
        g = Gases()
        eq = EquilibriumConstants(T = T)
    end
    @variables begin
        NH_eq(t), [unit = u"mol/m^3", description = "Equilibrium NH3 + NH4", guess = 1e-8]
        Na_eq(t), [unit = u"mol/m^3", description = "Equilibrium Na+", guess = 1e-8]
        Ca_eq(t), [unit = u"mol/m^3", description = "Equilibrium Ca2+", guess = 1e-8]
        K_eq(t), [unit = u"mol/m^3", description = "Equilibrium K+", guess = 1e-8]
        Mg_eq(t), [unit = u"mol/m^3", description = "Equilibrium Mg2+", guess = 1e-8]
        Cl_eq(t), [unit = u"mol/m^3", description = "Equilibrium Cl-", guess = 1e-8]
        NO3_eq(t), [unit = u"mol/m^3", description = "Equilibrium NO3-", guess = 1e-8]
        SO4_eq(t), [unit = u"mol/m^3", description = "Equilibrium SO4 2-", guess = 1e-8]

        TotalNH(t), [unit = u"mol/m^3", description = "Total NH3 + NH4", guess = 1e-8]
        TotalNa(t), [unit = u"mol/m^3", description = "Total Na+", guess = 1e-8]
        TotalCa(t), [unit = u"mol/m^3", description = "Total Ca2+", guess = 1e-8]
        TotalK(t), [unit = u"mol/m^3", description = "Total K+", guess = 1e-8]
        TotalMg(t), [unit = u"mol/m^3", description = "Total Mg2+", guess = 1e-8]
        TotalCl(t), [unit = u"mol/m^3", description = "Total Cl-", guess = 1e-8]
        TotalNO3(t), [unit = u"mol/m^3", description = "Total NO3-", guess = 1e-8]
        TotalSO4(t), [unit = u"mol/m^3", description = "Total SO4 2-", guess = 1e-8]
        mass_total(t), [unit = u"mol/m^3", description = "Total mass in the system"]

        NH_extra(t), [unit = u"mol/m^3", description = "Extra NH above equilibrium"]
        Na_extra(t), [unit = u"mol/m^3", description = "Extra Na above equilibrium"]
        Ca_extra(t), [unit = u"mol/m^3", description = "Extra Ca above equilibrium"]
        K_extra(t), [unit = u"mol/m^3", description = "Extra K above equilibrium"]
        Mg_extra(t), [unit = u"mol/m^3", description = "Extra Mg above equilibrium"]
        Cl_extra(t), [unit = u"mol/m^3", description = "Extra Cl above equilibrium"]
        NO3_extra(t), [unit = u"mol/m^3", description = "Extra NO3 above equilibrium"]
        SO4_extra(t), [unit = u"mol/m^3", description = "Extra SO4 above equilibrium"]

        R_1(t), [description = "Total sulfate ratio from (Section 3.1)"]
        R_2(t), [description = "Crustal species and sodium ratio (Section 3.1)"]
        R_3(t), [description = "Crustal species ratio (Section 3.1)"]

        type1(t), [description = "Sulfate rich (free acid) aerosol type (Table 3)"]
        type2(t), [description = "Sulfate rich aerosol type (Table 3)"]
        type3(t), [description = "Sulfate poor, crustal & sodium poor aerosol (Table 3)"]
        type4(t), [description = "Sulfate poor, crustal & sodium rich, crustal poor"]
        type5(t), [description = "Sulfate poor, crustal & sodium rich, crustal rich"]
    end
    @equations begin
        # Reactions based on information in Table 2 of Fountoukis and Nenes (2007).
        # The left-hand side of the reaction equation is the equilibrium constant and the
        # right-hand side is the ratio of the product and reactant activities.
        eq.r1.logK_eq ~ aq.CaNO32.loga
        eq.r2.logK_eq ~ aq.CaCl2.loga
        eq.r3.logK_eq ~ aq.CaSO4.loga + 2log(RH)
        eq.r4.logK_eq ~ aq.K2SO4.loga
        eq.r5.logK_eq ~ aq.KHSO4.loga
        eq.r6.logK_eq ~ aq.KNO3.loga
        eq.r7.logK_eq ~ aq.KCl.loga
        eq.r8.logK_eq ~ aq.MgSO4.loga
        eq.r9.logK_eq ~ aq.MgNO32.loga
        eq.r10.logK_eq ~ aq.MgCl2.loga
        eq.r11.logK_eq ~ 2log(aq.HSO4_dissociated.a / m_one) - log(aq.HSO4.a / m_one)
        eq.r12.logK_eq ~ log(aq.NH3.a / m_one) - log(g.NH3.p / p_one)
        eq.r13.logK_eq ~ 2log(aq.NH3_dissociated.a / m_one) -
            log(aq.NH3.a / m_one) - log(RH)
        eq.r14.logK_eq ~ aq.HNO3.loga - log(g.HNO3.p / p_one) # K1
        eq.r15.logK_eq ~ log(aq.HNO3_aq.a / m_one) - log(g.HNO3.p / p_one) # K1a
        #eq.r14.logK_eq - eq.r15.logK_eq ~ aq.HNO3.loga - log(aq.HNO3_aq.a / m_one) # K1b, from Table 2 footnote ♠
        eq.r16.logK_eq ~ aq.HCl.loga - log(g.HCl.p / p_one) # K2
        eq.r17.logK_eq ~ log(aq.HCl_aq.a / m_one) - log(g.HCl.p / p_one) # K2a
        #eq.r16.logK_eq - eq.r17.logK_eq ~ aq.HCl.loga - log(aq.HCl_aq.a / m_one) # K2b, from Table 2 footnote ♦
        eq.r18.logK_eq ~ 2log(aq.H2O_dissociated.a / m_one) - log(RH)
        eq.r19.logK_eq ~ aq.Na2SO4.loga
        eq.r20.logK_eq ~ aq.NH42SO4.loga
        eq.r21.logK_eq ~ log(g.NH3.p / p_one) + log(g.HCl.p / p_one)
        eq.r22.logK_eq ~ aq.NaNO3.loga
        eq.r23.logK_eq ~ aq.NaCl.loga
        eq.r24.logK_eq ~ aq.NaHSO4.loga
        eq.r25.logK_eq ~ log(g.NH3.p / p_one) + log(g.HNO3.p / p_one)
        eq.r26.logK_eq ~ aq.NH4HSO4.loga
        eq.r27.logK_eq ~ aq.NH43HSO42.loga

        g.HCl.M ~ g.NH3.M

        # Aqueous equilibrium mass Balance
        NH_eq ~ aq.NH4.M + aq.NH3.M #+ g.NH3.M
        Na_eq ~ aq.Na.M
        Ca_eq ~ aq.Ca.M
        K_eq ~ aq.K.M
        Mg_eq ~ aq.Mg.M
        Cl_eq ~ aq.Cl.M + aq.HCl_aq.M #+ g.HCl.M
        NO3_eq ~ aq.NO3.M + aq.HNO3_aq.M #+ g.HNO3.M
        SO4_eq ~ aq.SO4.M + aq.HSO4.M #+ aq.HSO4_dissociated.M

        NH_extra ~ TotalNH - NH_eq
        Na_extra ~ TotalNa - Na_eq
        Ca_extra ~ TotalCa - Ca_eq
        K_extra ~ TotalK - K_eq
        Mg_extra ~ TotalMg - Mg_eq
        Cl_extra ~ TotalCl - Cl_eq
        NO3_extra ~ TotalNO3 - NO3_eq
        SO4_extra ~ TotalSO4 - SO4_eq
        M_zero ~ min(NH_extra, Na_extra, Ca_extra, K_extra, Mg_extra,
            Cl_extra, NO3_extra, SO4_extra)

        D(TotalNH) ~ 0.0
        D(TotalNa) ~ 0.0 # 0.2*TotalNa/M_one * cos(0.5π * t/t_one) * M_one / t_one
        D(TotalCa) ~ 0.0
        D(TotalK) ~ 0.0
        D(TotalMg) ~ 0.0
        D(TotalCl) ~ 0.0
        D(TotalNO3) ~ 0.0 # 0.2*TotalNO3/M_one * cos(3π * t/t_one) * M_one / t_one
        D(TotalSO4) ~ 0.0 # 0.2*TotalSO4/M_one * sin(2π * t/t_one) * M_one / t_one

        # Aerosol types from Section 3.1 and Table 3
        R_1 ~ (TotalNH + TotalCa + TotalK + TotalMg + TotalNa) / TotalSO4
        R_2 ~ (TotalCa + TotalK + TotalMg + TotalNa) / TotalSO4
        R_3 ~ (TotalCa + TotalK + TotalMg) / TotalSO4
        type1 ~ 1 - (tanh((R_1 - 1) * 30) + 1) / 2
        type2 ~ min((tanh((R_1 - 1) * 30) + 1) / 2, 1-(tanh((R_1 - 2) * 30) + 1) / 2)
        type3 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, 1-(tanh((R_2 - 2) * 30) + 1) / 2)
        type4 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, (tanh((R_2 - 2) * 30) + 1) / 2,
            1-(tanh((R_3 - 2) * 30) + 1) / 2)
        type5 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, (tanh((R_2 - 2) * 30) + 1) / 2,
            1-(tanh((R_3 - 2) * 30) + 1) / 2)
    end
end

@doc """
    Isorropia(t)
    Isorropia(t, 1:3) # Only include the first three of the 27 reactions.
    Isorropia(t, 1) # Only include the first reaction.
    Isorropia(t, :all) # Choose a pre-specified set of reactions. Current options are `:all` and `:type1` (referring to first aerosol type in Fountoukis and Nenes (2007) Table 3).

An implementation of ISORROPIA II, a model for the thermodynamic equilibrium of gas-aerosol interactions, as described in:

> Fountoukis, C. and Nenes, A., 2007. ISORROPIA II: a computationally efficient thermodynamic equilibrium model for K+–Ca 2+–Mg 2+–NH 4+–Na+–SO 4 2−–NO 3−–Cl−–H 2 O aerosols. Atmospheric Chemistry and Physics, 7(17), pp.4639-4659.
""" Isorropia

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

# Ratios from Section 3.1
# R_1 ~ (totalNH + totalCa + totalK + totalMg + totalNa) / totalSO4
# R_2 ~ (totalCa + totalK + totalMg + totalNa) / totalSO4
# R_3 ~ (totalCa + totalK + totalMg) / totalSO4

spc = [:CaNO32s
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
       :NH43HSO42aq]

unique(spc)

end
