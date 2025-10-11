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
    @components begin
        aq = Aqueous(T = T, RH = RH)
       # s = Solids()
       #g = Gases()
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
        # eq.r1.logK_eq ~ aq.CaNO32.loga
        eq.r2.logK_eq ~ aq.CaCl2.loga
        # eq.r3.logK_eq ~ aq.CaSO4.loga * RH^2 # Not using because CaSO4 precipitates completely (Table 4 footnote a).
        eq.r4.logK_eq ~ aq.K2SO4.loga
        #eq.r5.logK_eq ~ aq.KHSO4.loga
        # eq.r6.logK_eq ~ aq.KNO3.loga
        # eq.r7.logK_eq ~ aq.KCl.loga
        eq.r8.logK_eq ~ aq.MgSO4.loga
        eq.r9.logK_eq ~ aq.MgNO32.loga
        eq.r10.logK_eq ~ aq.MgCl2.loga
        #eq.r11.logK_eq ~ log(aq.H.a) + log(aq.SO4.a) - log(aq.HSO4.a)
     #   eq.r12.logK_eq ~ log(aq.NH3.a / m_one) - log(g.NH3.p / p_one)
    #    eq.r13.logK_eq ~ log(aq.NH4.a) + log(aq.OH.a) - log(aq.NH3.a) - log(RH)
        #eq.r14.logK_eq ~ aq.HNO3.loga - log(g.HNO3.p) # K1
        #eq.r15.logK_eq ~ log(aq.HNO3_aq.a) - log(g.HNO3.p) # K1a
        eq.r14.logK_eq - eq.r15.logK_eq ~ aq.HNO3.loga - log(aq.HNO3_aq.a / m_one) # K1b, from Table 2 footnote ♠
        # eq.r16.logK_eq ~ aq.HCl.loga - log(g.HCl.p) # K2
        # eq.r17.logK_eq ~ log(aq.HCl_aq.a) - log( g.HCl.p) # K2a
        #eq.r16.logK_eq - eq.r17.logK_eq ~ aq.HCl.loga - log(aq.HCl_aq.a) # K2b, from Table 2 footnote ♦
        #eq.r18.logK_eq ~ log(aq.H.a) + log(aq.OH.a) - log(RH)
        eq.r19.logK_eq ~ aq.Na2SO4.loga
        eq.r20.logK_eq ~ aq.NH42SO4.loga
    #    eq.r21.logK_eq ~ log(g.NH3.p / p_one) + log(g.HCl.p / p_one)
        #eq.r22.logK_eq ~ aq.NaNO3.loga
        #eq.r23.logK_eq ~ aq.NaCl.loga
        #eq.r24.logK_eq ~ aq.NaHSO4.loga
        #eq.r25.logK_eq ~ log(g.NH3.p) + log(g.HNO3.p)
        eq.r26.logK_eq ~ aq.NH4HSO4.loga
        eq.r27.logK_eq ~ aq.NH43HSO42.loga

        # Mass Balance
        # TotalNH ~ (aq.NH4.m  + aq.NH3.m) * aq.W #+ g.NH3.M #+ s.NH4 #
        # TotalNa ~ aq.Na.m * aq.W #+ s.Na
        # TotalCa ~ aq.Ca.m * aq.W #+ s.Ca
        # TotalK ~ aq.K.m * aq.W #+ s.K
        # TotalMg ~ aq.Mg.m * aq.W #+ s.Mg
        # TotalCl ~  aq.Cl.m * aq.W #+ g.HCl.M #+ s.Cl #
        # TotalNO3 ~ aq.NO3.m * aq.W #+ g.HNO3.M #+ s.NO3 #
        # TotalSO4 ~ (aq.SO4.m + aq.HSO4.m) * aq.W #+ g.H2SO4.M #+ s.SO4 + s.HSO4 #

        aq.NH3.m ~ m_one
        aq.HCl_aq.m ~ m_one
        aq.HNO3_aq.m ~ m_one

        # g.NH3.M ~ 0.0
        # g.HNO3.M ~ 0.0
        # g.HCl.M ~ M_one

        D(aq.NH42SO4.M) ~ 0.0
        #D(TotalNH) ~ 0.0
        # D(TotalNa) ~ 0.0
        # D(TotalCa) ~ 0.0
        # D(TotalK) ~ 0.0
        # D(TotalMg) ~ 0.0
        # D(TotalCl) ~ 0.0
        # D(TotalNO3) ~ 0.0
        # D(TotalSO4) ~ 0.0

        # aq.NH4Cl.M ~ 0 # Not present in Table 2?
        # aq.NH4NO3.M ~ 0 # Not present in Table 2?

        # aq.Mg.m ~ 0
        # aq.OH.m ~ 0
        # aq.CaNO32.M ~ 0
        # aq.CaCl2.M ~ 0
        # aq.K2SO4.M ~ 0
        # aq.KNO3.M ~ 0
        # aq.KCl.M ~ 0
        # aq.MgSO4.M ~ 0
        # aq.MgNO32.M ~ 0
        # aq.MgCl2.M ~ 0
        # #aq.NaCl.M ~ 0
        # # aq.Na2SO4.M ~ 0
        # aq.NaNO3.M ~ 0
        # aq.NH42SO4.M ~ 0
        # aq.NH4NO3.M ~ 0
        # aq.NH4Cl.M ~ 0
        # aq.NH43HSO42.M ~ 0
        # # aq.HNO3.M ~ 0
        # #aq.HCl.M ~ 0
        # s.CaNO32.M ~ 0
        # s.CaCl2.M ~ 0
        # s.K2SO4.M ~ 0
        # s.KNO3.M ~ 0
        # s.KCl.M ~ 0
        # #s.MgSO4.M ~ 0
        # s.MgNO32.M ~ 0
        # s.MgCl2.M ~ 0
        # s.NaCl.M ~ 0
        # s.NaNO3.M ~ 0
        # s.Na2SO4.M ~ 0
        # s.NH4Cl.M ~ 0
        # s.NH4NO3.M ~ 0
        # #s.NH42SO4.M ~ 0
        # s.NH43HSO42.M ~ 0
        # g.HNO3.M ~ 0
        # g.HCl.M ~ 0
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
