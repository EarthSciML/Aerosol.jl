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
        NH_eq(t), [unit = u"mol/kg", description = "Equilibrium NH3 + NH4", guess = 1]
        Na_eq(t), [unit = u"mol/kg", description = "Equilibrium Na+", guess = 1]
        Ca_eq(t), [unit = u"mol/kg", description = "Equilibrium Ca2+", guess = 1]
        K_eq(t), [unit = u"mol/kg", description = "Equilibrium K+", guess = 1]
        Mg_eq(t), [unit = u"mol/kg", description = "Equilibrium Mg2+", guess = 1]
        Cl_eq(t), [unit = u"mol/kg", description = "Equilibrium Cl-", guess = 1]
        NO3_eq(t), [unit = u"mol/kg", description = "Equilibrium NO3-", guess = 1]
        SO4_eq(t), [unit = u"mol/kg", description = "Equilibrium SO4 2-", guess = 1]

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
        eq.r1.logK_eq ~ aq.CaNO32.loga_eq
        eq.r2.logK_eq ~ aq.CaCl2.loga_eq
        eq.r3.logK_eq ~ aq.CaSO4.loga_eq + 2log(RH)
        eq.r4.logK_eq ~ aq.K2SO4.loga_eq
        eq.r5.logK_eq ~ aq.KHSO4.loga_eq
        eq.r6.logK_eq ~ aq.KNO3.loga_eq
        eq.r7.logK_eq ~ aq.KCl.loga_eq
        eq.r8.logK_eq ~ aq.MgSO4.loga_eq
        eq.r9.logK_eq ~ aq.MgNO32.loga_eq
        eq.r10.logK_eq ~ aq.MgCl2.loga_eq
        eq.r12.logK_eq ~ log(aq.NH3.m_aq / m_one) - log(g.NH3.p / p_one) - 13 # FIXME(CT): Added -13 to get reasonable results; not sure why this is needed.
        eq.r13.logK_eq ~ aq.NH3_dissociated.loga_aq - log(aq.NH3.m_aq / m_one) + log(RH) - 8 # FIXME(CT): Added -8 to get reasonable results; not sure why this is needed.
        eq.r14.logK_eq ~ aq.HNO3.loga_aq - log(g.HNO3.p / p_one) - 16 # FIXME(CT): Added -16 to get reasonable results; not sure why this is needed.
        eq.r15.logK_eq ~ log(aq.HNO3_aq.m_aq / m_one) - log(g.HNO3.p / p_one)
        eq.r16.logK_eq ~ aq.HCl.loga_aq - log(g.HCl.p / p_one) - 23 # FIXME(CT): Added -23 to get reasonable results; not sure why this is needed.
        eq.r17.logK_eq ~ log(aq.HCl_aq.m_aq / m_one) - log(g.HCl.p / p_one)
        eq.r18.logK_eq ~ aq.H2O_dissociated.loga_aq - log(RH)
        eq.r19.logK_eq ~ aq.Na2SO4.loga_eq
        eq.r20.logK_eq ~ aq.NH42SO4.loga_eq
        eq.r22.logK_eq ~ aq.NaNO3.loga_eq
        eq.r23.logK_eq ~ aq.NaCl.loga_eq
        eq.r24.logK_eq ~ aq.NaHSO4.loga_eq
        eq.r26.logK_eq ~ aq.NH4HSO4.loga_eq
        eq.r27.logK_eq ~ aq.NH43HSO42.loga_eq

        #eq.r11.logK_eq ~ 2log(aq.HSO4_dissociated.a / m_one) - aq.H2SO4.loga # H2SO4 fully dissociates, so the concentration of HSO4 is the same as the concentration of HHSO4.
        #eq.r21.logK_eq ~ log(g.NH3.p / p_one) + log(g.HCl.p / p_one)
        # eq.r25.logK_eq ~ log(g.NH3.p / p_one) + log(g.HNO3.p / p_one)

        #aq.NH3.m ~ 0
       #  aq.HCl_aq.m_eq ~ 1e-20 * m_one
    #     aq.HNO3_aq.m_eq ~ 1e-20 * m_one
    #     aq.H2O_dissociated.m_eq ~ 1e-20 * m_one
 #         g.HNO3.M ~ 0
 # g.HCl.M ~ 0
 #g.NH3.M ~ 0

        # Aqueous equilibrium mass Balance
        NH_eq ~ aq.NH3_dissociated.m_eq + aq.NH3.m_eq #+ g.NH3.M ## XXXXXXXX
        Na_eq ~ aq.Na.m_eq
        Ca_eq ~ aq.Ca.m_eq
        K_eq ~ aq.K.m_eq
        Mg_eq ~ aq.Mg.m_eq
        Cl_eq ~ aq.HCl.m_eq + aq.HCl_aq.m_eq #+ g.HCl.M
        NO3_eq ~ aq.HNO3.m_eq + aq.HNO3_aq.m_eq #+ g.HNO3.M
        SO4_eq ~ aq.HSO4_dissociated.m_eq + aq.HHSO4.m_eq

#        D(Na_eq) ~ 0.0
        NH_extra ~ TotalNH - g.NH3.M - NH_eq * aq.W_eq
        Na_extra ~ TotalNa - Na_eq * aq.W_eq
        Ca_extra ~ TotalCa - Ca_eq * aq.W_eq
        K_extra ~ TotalK - K_eq * aq.W_eq
        Mg_extra ~ TotalMg - Mg_eq * aq.W_eq
        Cl_extra ~ TotalCl - g.HCl.M - Cl_eq * aq.W_eq
        NO3_extra ~ TotalNO3 - g.HNO3.M - NO3_eq * aq.W_eq
        SO4_extra ~ TotalSO4 - SO4_eq * aq.W_eq
        aq.W_eq ~ min(
            (TotalNH - g.NH3.M) / NH_eq,
            TotalNa / Na_eq,
            TotalCa / Ca_eq,
            TotalK / K_eq,
            TotalMg / Mg_eq,
            (TotalCl - g.HCl.M) / Cl_eq,
            (TotalNO3 - g.HNO3.M) / NO3_eq,
            TotalSO4 / SO4_eq,
        )
        # M_zero ~ min(NH_extra, Na_extra, Ca_extra, K_extra, Mg_extra,
        #     Cl_extra, NO3_extra, SO4_extra)

        D(TotalNH) ~ 0.0
        D(TotalNa) ~ 0.0 # 0.2*TotalNa/M_one * cos(0.5π * t/t_one) * M_one / t_one
        D(TotalCa) ~ 0.0
        D(TotalK) ~ 0.0
        D(TotalMg) ~ 0.0
        D(TotalCl) ~ 0.0
        D(TotalNO3) ~ 0.0 # 0.2*TotalNO3/M_one * cos(3π * t/t_one) * M_one / t_one
        D(TotalSO4) ~ 0.0 # 0.2*TotalSO4/M_one * sin(2π * t/t_one) * M_one / t_one

        # # Aerosol types from Section 3.1 and Table 3
        # R_1 ~ (TotalNH + TotalCa + TotalK + TotalMg + TotalNa) / TotalSO4
        # R_2 ~ (TotalCa + TotalK + TotalMg + TotalNa) / TotalSO4
        # R_3 ~ (TotalCa + TotalK + TotalMg) / TotalSO4
        # type1 ~ 1 - (tanh((R_1 - 1) * 30) + 1) / 2
        # type2 ~ min((tanh((R_1 - 1) * 30) + 1) / 2, 1-(tanh((R_1 - 2) * 30) + 1) / 2)
        # type3 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, 1-(tanh((R_2 - 2) * 30) + 1) / 2)
        # type4 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, (tanh((R_2 - 2) * 30) + 1) / 2,
        #     1-(tanh((R_3 - 2) * 30) + 1) / 2)
        # type5 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, (tanh((R_2 - 2) * 30) + 1) / 2,
        #     1-(tanh((R_3 - 2) * 30) + 1) / 2)
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
