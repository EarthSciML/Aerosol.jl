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
    end
    @components begin
        aq = Aqueous(T = T, RH = RH)
        s = Solids()
        g = Gases()
        eq = EquilibriumConstants(T = T)
    end
    @variables begin
        aqNH(t), [unit = u"mol/m^3", description = "Aqueous NH3 + NH4+", guess = 1e-20]
        sNH(t), [unit = u"mol/m^3", description = "Solid NH3 + NH4+", guess = 1e-20]
        aqNa(t), [unit = u"mol/m^3", description = "Aqueous Na+", guess = 1e-20]
        sNa(t), [unit = u"mol/m^3", description = "Solid Na+", guess = 1e-20]
        aqCa(t), [unit = u"mol/m^3", description = "Aqueous Ca2+", guess = 1e-20]
        sCa(t), [unit = u"mol/m^3", description = "Solid Ca2+", guess = 1e-20]
        aqK(t), [unit = u"mol/m^3", description = "Aqueous K+", guess = 1e-20]
        sK(t), [unit = u"mol/m^3", description = "Solid K+", guess = 1e-20]
        aqMg(t), [unit = u"mol/m^3", description = "Aqueous Mg2+", guess = 1e-20]
        sMg(t), [unit = u"mol/m^3", description = "Solid Mg2+", guess = 1e-20]
        aqCl(t), [unit = u"mol/m^3", description = "Aqueous Cl-", guess = 1e-20]
        sCl(t), [unit = u"mol/m^3", description = "Solid Cl-", guess = 1e-20]
        aqNO3(t), [unit = u"mol/m^3", description = "Aqueous NO3-", guess = 1e-20]
        sNO3(t), [unit = u"mol/m^3", description = "Solid NO3-", guess = 1e-20]
        aqSO4(t), [unit = u"mol/m^3", description = "Aqueous SO4 2-", guess = 1e-20]
        sSO4(t), [unit = u"mol/m^3", description = "Solid SO4 2-", guess = 1e-20]

        #! format: off
        TotalNH(t), [unit = u"mol/m^3", description = "Total NH3 + NH4 Molarity", guess = 1e-20]
        TotalNa(t), [unit = u"mol/m^3", description = "Total Na+ Molarity", guess = 1e-20]
        TotalCa(t), [unit = u"mol/m^3", description = "Total Ca2+ Molarity", guess = 1e-20]
        TotalK(t), [unit = u"mol/m^3", description = "Total K+ Molarity", guess = 1e-20]
        TotalMg(t), [unit = u"mol/m^3", description = "Total Mg2+ Molarity", guess = 1e-20]
        TotalCl(t), [unit = u"mol/m^3", description = "Total Cl- Molarity", guess = 1e-20]
        TotalNO3(t), [unit = u"mol/m^3", description = "Total NO3- Molarity", guess = 1e-20]
        TotalSO4(t), [unit = u"mol/m^3", description = "Total SO4 2- Molarity", guess = 1e-20]
        #! format: on
    end
    @equations begin
        # Reactions based on information in Table 2 of Fountoukis and Nenes (2007).
        # The left-hand side of the reaction equation is the equilibrium constant and the
        # right-hand side is the ratio of the product and reactant activities.
        # eq.r1.K_eq * eq.k1_unit ~ aq.a_CaNO32
        # eq.r2.K_eq * eq.k2_unit ~ aq.a_CaCl2
        # eq.r3.K_eq * eq.k3_unit ~ aq.a_CaSO4 * RH^2 # Not using in favor of equation below
        aq.CaSO4.M ~ 0.0 # From Table 4 footnote a, CaSO4 has an activity coefficient of zero and therefore 0 concentration.
        # eq.r4.K_eq * eq.k4_unit ~ aq.a_K2SO4
        #eq.r5.K_eq * eq.k5_unit ~ aq.a_KHSO4
        # eq.r6.K_eq * eq.k6_unit ~ aq.a_KNO3
        # eq.r7.K_eq * eq.k7_unit ~ aq.a_KCl
        # eq.r8.K_eq * eq.k8_unit ~ aq.a_MgSO4
        # eq.r9.K_eq * eq.k9_unit ~ aq.a_MgNO32
        # eq.r10.K_eq * eq.k10_unit ~ aq.a_MgCl2
        eq.r11.K_eq * eq.k11_unit ~ aq.H.a * aq.SO4.a / aq.HSO4.a
        eq.r12.K_eq * eq.k12_unit ~ aq.NH3.a / g.NH3.p
        eq.r13.K_eq * eq.k13_unit ~ aq.NH4.a * aq.OH.a / aq.NH3.a / RH
        # eq.r14.K_eq * eq.k14_unit ~ aq.a_HNO3 / g.HNO3.p # K1
        # eq.r15.K_eq * eq.k15_unit ~ aq.HNO3_aq.a / g.HNO3.p # K1a
        eq.r14.K_eq / eq.r15.K_eq * eq.k1b_unit ~ aq.a_HNO3 / aq.HNO3_aq.a # K1b, from Table 2 footnote ♠
        # eq.r16.K_eq * eq.k16_unit ~ aq.a_HCl / g.HCl.p # K2
        # eq.r17.K_eq * eq.k17_unit ~ aq.HCl_aq.a / g.HCl.p # K2a
        eq.r16.K_eq / eq.r17.K_eq * eq.k2b_unit ~ aq.a_HCl / aq.HCl_aq.a # K2b, from Table 2 footnote ♦
        eq.r18.K_eq * eq.k18_unit ~ aq.H.a * aq.OH.a / RH
        # eq.r19.K_eq * eq.k19_unit ~ aq.a_Na2SO4
        # eq.r20.K_eq * eq.k20_unit ~ aq.a_NH42SO4
        # eq.r21.K_eq * eq.k21_unit ~ g.NH3.p * g.HCl.p
        #eq.r22.K_eq * eq.k22_unit ~ aq.a_NaNO3
        #eq.r23.K_eq * eq.k23_unit ~ aq.a_NaCl
        eq.r24.K_eq * eq.k24_unit ~ aq.a_NaHSO4
        # eq.r25.K_eq * eq.k25_unit ~ g.NH3.p * g.HNO3.p
        eq.r26.K_eq * eq.k26_unit ~ aq.a_NH4HSO4
        #eq.r27.K_eq * eq.k27_unit ~ aq.a_NH43HSO42

        # Mass Balance
        aqNH ~ (aq.NH3.m + aq.NH4.m) * aq.W
        sNH ~ s.NH4Cl.M + s.NH4NO3.M + 2s.NH42SO4.M + s.NH4HSO4.M + 3s.NH43HSO42.M
        TotalNH ~ g.NH3.M + aqNH + sNH

        aqNa ~ aq.Na.m * aq.W
        sNa ~ s.NaCl.M + s.NaNO3.M + 2s.Na2SO4.M + s.NaHSO4.M
        TotalNa ~ aqNa + sNa

        aqCa ~ aq.Ca.m * aq.W
        sCa ~ s.CaNO32.M + s.CaCl2.M + s.CaSO4.M
        TotalCa ~ aqCa + sCa

        aqK ~ aq.K.m * aq.W
        sK ~ s.KHSO4.M + 2s.K2SO4.M + s.KNO3.M + s.KCl.M
        TotalK ~ aqK + sK

        aqMg ~ aq.Mg.m * aq.W
        sMg ~ s.MgSO4.M + s.MgNO32.M + s.MgCl2.M
        TotalMg ~ aqMg + sMg

        aqCl ~ (aq.HCl_aq.m + aq.Cl.m) * aq.W
        sCl ~ 2s.CaCl2.M + s.KCl.M + 2s.MgCl2.M + s.NaCl.M + s.NH4Cl.M
        TotalCl ~ g.HCl.M + aqCl + sCl

        aqNO3 ~ (aq.HNO3_aq.m + aq.NO3.m) * aq.W
        sNO3 ~ 2s.CaNO32.M + s.KNO3.M + 2s.MgNO32.M + s.NaNO3.M + s.NH4NO3.M
        TotalNO3 ~ g.HNO3.M + aqNO3 + sNO3

        aqSO4 ~ (aq.SO4.m + aq.HSO4.m) * aq.W
        sSO4 ~ s.CaSO4.M + s.KHSO4.M + s.K2SO4.M + s.MgSO4.M + s.Na2SO4.M + s.NaHSO4.M +
               s.NH42SO4.M + s.NH4HSO4.M + 2s.NH43HSO42.M
        TotalSO4 ~ g.H2SO4.M + aqSO4 + sSO4

        D(TotalNH) ~ 0.0
        D(TotalNa) ~ 0.0
        D(TotalCa) ~ 0.0
        D(TotalK) ~ 0.0
        D(TotalMg) ~ 0.0
        D(TotalCl) ~ 0.0
        D(TotalNO3) ~ 0.0
        D(TotalSO4) ~ 0.0

        aq.NH4Cl.M ~ 0 # Not present in Table 2?
        aq.NH4NO3.M ~ 0 # Not present in Table 2?

        aq.Mg.m ~ 0
        aq.OH.m ~ 0
        aq.CaNO32.M ~ 0
        aq.CaCl2.M ~ 0
        aq.K2SO4.M ~ 0
        aq.KNO3.M ~ 0
        aq.KCl.M ~ 0
        aq.MgSO4.M ~ 0
        aq.MgNO32.M ~ 0
        aq.MgCl2.M ~ 0
        #aq.NaCl.M ~ 0
        # aq.Na2SO4.M ~ 0
        aq.NaNO3.M ~ 0
        aq.NH42SO4.M ~ 0
        aq.NH4NO3.M ~ 0
        aq.NH4Cl.M ~ 0
        aq.NH43HSO42.M ~ 0
        # aq.HNO3.M ~ 0
        #aq.HCl.M ~ 0
        s.CaNO32.M ~ 0
        s.CaCl2.M ~ 0
        s.K2SO4.M ~ 0
        s.KNO3.M ~ 0
        s.KCl.M ~ 0
        #s.MgSO4.M ~ 0
        s.MgNO32.M ~ 0
        s.MgCl2.M ~ 0
        s.NaCl.M ~ 0
        s.NaNO3.M ~ 0
        s.Na2SO4.M ~ 0
        s.NH4Cl.M ~ 0
        s.NH4NO3.M ~ 0
        #s.NH42SO4.M ~ 0
        s.NH43HSO42.M ~ 0
        g.HNO3.M ~ 0
        g.HCl.M ~ 0
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
