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
    @constants begin
        m_one = 1.0, [unit = u"mol/kg"]
        M_one = 1.0, [unit = u"mol/m^3"]
        p_one = 1.0, [unit = u"Constants.atm"]
        M_zero = 0.0, [unit = u"mol/m^3"]
        m_zero = 0.0, [unit = u"mol/kg"]
        w_inf = 1e30, [unit = u"kg/m^3"]
    end
    @components begin
        aq = Aqueous(T = T, RH = RH)
        #s = Solids()
        g = Gases()
        eq = EquilibriumConstants(T = T)
    end
    @variables begin
        NH_eq(t), [unit = u"mol/kg", description = "Equilibrium NH3 + NH4", guess = 19]
        Na_eq(t), [unit = u"mol/kg", description = "Equilibrium Na+", guess = 159]
        Ca_eq(t), [unit = u"mol/kg", description = "Equilibrium Ca2+", guess = 12]
        K_eq(t), [unit = u"mol/kg", description = "Equilibrium K+", guess = 6]
        Mg_eq(t), [unit = u"mol/kg", description = "Equilibrium Mg2+", guess = 134]
        Cl_eq(t), [unit = u"mol/kg", description = "Equilibrium Cl-", guess = 177]
        NO3_eq(t), [unit = u"mol/kg", description = "Equilibrium NO3-", guess = 116]
        SO4_eq(t), [unit = u"mol/kg", description = "Equilibrium SO4 2-", guess = 7]

        TotalNH(t)=2e-7, [unit = u"mol/m^3", description = "Total NH3 + NH4"]
        TotalNa(t)=4e-10, [unit = u"mol/m^3", description = "Total Na+"]
        TotalCa(t)=1e-8, [unit = u"mol/m^3", description = "Total Ca2+"]
        TotalK(t)=8.4e-9, [unit = u"mol/m^3", description = "Total K+"]
        TotalMg(t)=4.1e-10, [unit = u"mol/m^3", description = "Total Mg2+"]
        TotalCl(t)=2.7e-10, [unit = u"mol/m^3", description = "Total Cl-"]
        TotalNO3(t)=3.2e-8, [unit = u"mol/m^3", description = "Total NO3-"]
        TotalSO4(t)=1e-7, [unit = u"mol/m^3", description = "Total SO4 2-"]
        mass_total(t), [unit = u"mol/m^3", description = "Total mass in the system"]

        #! format: off
        NH_extra(t), [unit = u"mol/m^3", description = "Extra NH above equilibrium for allocation to salts", guess=1.4e-7]
        Na_extra(t), [unit = u"mol/m^3", description = "Extra Na above equilibrium for allocation to salts", guess=1.4e-7]
        Ca_extra(t), [unit = u"mol/m^3", description = "Extra Ca above equilibrium for allocation to salts"]
        K_extra(t), [unit = u"mol/m^3", description = "Extra K above equilibrium for allocation to salts"]
        Mg_extra(t), [unit = u"mol/m^3", description = "Extra Mg above equilibrium for allocation to salts"]
        Cl_extra(t), [unit = u"mol/m^3", description = "Extra Cl above equilibrium for allocation to salts"]
        NO3_extra(t), [unit = u"mol/m^3", description = "Extra NO3 above equilibrium for allocation to salts"]
        SO4_extra(t), [unit = u"mol/m^3", description = "Extra SO4 above equilibrium for allocation to salts"]
        #! format: on

        SO4_resid(t), [unit = u"mol/m^3", description = "Residual SO4 after salt formation"]
        NO3_resid(t), [unit = u"mol/m^3", description = "Residual NO3 after salt formation"]
        NH_resid(t), [unit = u"mol/m^3", description = "Residual NH after salt formation"]
        Cl_resid(t), [unit = u"mol/m^3", description = "Residual Cl after salt formation"]

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
        eq.r12.logK_eq ~ log(abs(aq.NH3.m_aq) / m_one) - log(g.NH3.p / p_one) - 13 # FIXME(CT): Added -13 to get reasonable results; not sure why this is needed.
        eq.r13.logK_eq ~ aq.NH3_dissociated.loga_aq - log(abs(aq.NH3.m_aq) / m_one) + log(RH) - 8 # FIXME(CT): Added -8 to get reasonable results; not sure why this is needed.
        eq.r14.logK_eq ~ aq.HNO3.loga_aq - log(g.HNO3.p / p_one) - 18 # FIXME(CT): Added -18 to get reasonable results; not sure why this is needed.
        eq.r15.logK_eq ~ log(abs(aq.HNO3_aq.m_aq) / m_one) - log(g.HNO3.p / p_one)
        eq.r16.logK_eq ~ aq.HCl.loga_aq - log(g.HCl.p / p_one) - 23 # FIXME(CT): Added -23 to get reasonable results; not sure why this is needed.
        eq.r17.logK_eq ~ log(abs(aq.HCl_aq.m_aq) / m_one) - log(g.HCl.p / p_one)
        eq.r18.logK_eq ~ aq.H2O_dissociated.loga_aq - log(RH)
        eq.r19.logK_eq ~ aq.Na2SO4.loga_eq
        eq.r20.logK_eq ~ aq.NH42SO4.loga_eq
        eq.r22.logK_eq ~ aq.NaNO3.loga_eq
        eq.r23.logK_eq ~ aq.NaCl.loga_eq
        eq.r24.logK_eq ~ aq.NaHSO4.loga_eq
        eq.r26.logK_eq ~ aq.NH4HSO4.loga_eq
        eq.r27.logK_eq ~ aq.NH43HSO42.loga_eq

        # These equations would result in an over-specified equilibrium.
        # eq.r11.logK_eq ~ 2log(aq.HSO4_dissociated.a / m_one) - aq.H2SO4.loga # H2SO4 fully dissociates, so the concentration of HSO4 is the same as the concentration of HHSO4.
        # eq.r21.logK_eq ~ log(g.NH3.p / p_one) + log(g.HCl.p / p_one)
        # eq.r25.logK_eq ~ log(g.NH3.p / p_one) + log(g.HNO3.p / p_one)

        ## Aqueous equilibrium mass Balance
        # The gaseous concentrations should never be larger than the total concentrations,
        # but it is numerically possible so we handle that possibility to avoid negative
        # numbers.
        # NH_eq ~ ifelse(TotalNH > g.NH3.M, aq.NH3_dissociated.m_eq + aq.NH3.m_eq, m_zero)
        NH_eq ~ aq.NH3_dissociated.m_eq + aq.NH3.m_eq
        Na_eq ~ aq.Na.m_eq
        Ca_eq ~ aq.Ca.m_eq
        K_eq ~ aq.K.m_eq
        Mg_eq ~ aq.Mg.m_eq
        #Cl_eq ~ ifelse(TotalCl > g.HCl.M, aq.HCl.m_eq + aq.HCl_aq.m_eq, m_zero)
        Cl_eq ~ aq.HCl.m_eq + aq.HCl_aq.m_eq
        #NO3_eq ~ ifelse(TotalNO3 > g.HNO3.M, aq.HNO3.m_eq + aq.HNO3_aq.m_eq, m_zero)
        NO3_eq ~ aq.HNO3.m_eq + aq.HNO3_aq.m_eq
        SO4_eq ~ aq.HSO4_dissociated.m_eq + aq.HHSO4.m_eq

        NH_extra ~ TotalNH - g.NH3.M - NH_eq * aq.W_eq
        Na_extra ~ TotalNa - Na_eq * aq.W_eq
        Ca_extra ~ TotalCa - Ca_eq * aq.W_eq
        K_extra ~ TotalK - K_eq * aq.W_eq
        Mg_extra ~ TotalMg - Mg_eq * aq.W_eq
        Cl_extra ~ TotalCl - g.HCl.M - Cl_eq * aq.W_eq
        NO3_extra ~ TotalNO3 - g.HNO3.M - NO3_eq * aq.W_eq
        SO4_extra ~ TotalSO4 - SO4_eq * aq.W_eq
        aq.W_eq ~ min(
            #ifelse(TotalNH > g.NH3.M, (TotalNH - g.NH3.M) / NH_eq, w_inf),
            (TotalNH - g.NH3.M) / NH_eq,
            TotalNa / Na_eq,
            TotalCa / Ca_eq,
            TotalK / K_eq,
            TotalMg / Mg_eq,
            #ifelse(TotalCl > g.HCl.M, (TotalCl - g.HCl.M) / Cl_eq, w_inf),
            (TotalCl - g.HCl.M) / Cl_eq,
            #ifelse(TotalNO3 > g.HNO3.M, (TotalNO3 - g.HNO3.M) / NO3_eq, w_inf),
            (TotalNO3 - g.HNO3.M) / NO3_eq,
            TotalSO4 / SO4_eq
        )

        D(TotalNH) ~ 0.0
        D(TotalNa) ~ 0.0
        D(TotalCa) ~ 0.0
        D(TotalK) ~ 0.0
        D(TotalMg) ~ 0.0
        D(TotalCl) ~ 0.0
        D(TotalNO3) ~ 0.0
        D(TotalSO4) ~ 0.0

        # Aerosol types from Section 3.1 and Table 3
        R_1 ~ (TotalNH + TotalCa + TotalK + TotalMg + TotalNa) / TotalSO4
        R_2 ~ (TotalCa + TotalK + TotalMg + TotalNa) / TotalSO4
        R_3 ~ (TotalCa + TotalK + TotalMg) / TotalSO4
        type1 ~ 1 - (tanh((R_1 - 1) * 30) + 1) / 2
        type2 ~ min((tanh((R_1 - 1) * 30) + 1) / 2, 1 - (tanh((R_1 - 2) * 30) + 1) / 2)
        type3 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, 1 - (tanh((R_2 - 2) * 30) + 1) / 2)
        type4 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, (tanh((R_2 - 2) * 30) + 1) / 2,
            1 - (tanh((R_3 - 2) * 30) + 1) / 2)
        type5 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, (tanh((R_2 - 2) * 30) + 1) / 2,
            1 - (tanh((R_3 - 2) * 30) + 1) / 2)

        # Solid mass balances
        # First fill in salts that are needed to guarantee space, or that are otherwise prioritized.
        aq.Na2SO4.M_salt ~ min(Na_extra/2, SO4_extra) # Na preferentially forms Na2SO4 over other salts (section 3.2)
        aq.K2SO4.M_salt ~ min(K_extra/2, max(SO4_extra - aq.Na2SO4.M_salt, M_zero)) # K preferentially forms K2SO4 over other salts (section 3.2)
        aq.CaCl2.M_salt ~ min(Ca_extra, Cl_extra/2)
         aq.MgCl2.M_salt ~ min(Mg_extra, Cl_extra/2 - aq.CaCl2.M_salt)
         aq.MgNO32.M_salt ~ min(Mg_extra - aq.MgCl2.M_salt, NO3_extra/2)
        aq.NH4HSO4.M_salt ~ min(NH_extra, max(SO4_extra - aq.Na2SO4.M_salt -
            aq.K2SO4.M_salt, M_zero))
         aq.CaNO32.M_salt ~ min(Ca_extra, max(NO3_extra/2 - aq.CaNO32.M_salt, M_zero))
        aq.NaHSO4.M_salt ~ min(max(Na_extra - 2aq.Na2SO4.M_salt, M_zero),
            max(SO4_extra - aq.NH4HSO4.M_salt - aq.Na2SO4.M_salt - aq.K2SO4.M_salt, M_zero))
        # Next, fill in salts in order of increasing DRH unless otherwise noted.
        aq.NH4NO3.M_salt ~ min(max(NH_extra - aq.NH4HSO4.M_salt, M_zero),
            max(NO3_extra - 2aq.CaNO32.M_salt - 2aq.MgNO32.M_salt, M_zero))
        aq.NH43HSO42.M_salt ~ min(
            max((NH_extra  - aq.NH4HSO4.M_salt- aq.NH4NO3.M_salt )/3, M_zero),
            max((SO4_extra - aq.Na2SO4.M_salt- aq.K2SO4.M_salt - aq.NH4HSO4.M_salt -
                aq.NH4HSO4.M_salt)/2, M_zero))
        aq.NaNO3.M_salt ~ min(max(Na_extra - 2aq.Na2SO4.M_salt - aq.NaHSO4.M_salt, M_zero),
            max(NO3_extra - 2aq.CaNO32.M_salt - 2aq.MgNO32.M_salt - aq.NH4NO3.M_salt, M_zero))
        aq.NaCl.M_salt ~ min(max(Na_extra - 2aq.Na2SO4.M_salt - aq.NaHSO4.M_salt -
            aq.NaNO3.M_salt, M_zero),
            max(Cl_extra - 2aq.CaCl2.M_salt - 2aq.MgCl2.M_salt, M_zero))
        aq.NH4Cl.M_salt ~ min(NH_extra - aq.NH4HSO4.M_salt - aq.NH4NO3.M_salt -
            3aq.NH43HSO42.M_salt,
            Cl_extra - 2aq.CaCl2.M_salt - 2aq.MgCl2.M_salt - aq.NaCl.M_salt)
        aq.NH42SO4.M_salt ~ min(max((NH_extra - aq.NH4HSO4.M_salt - aq.NH4NO3.M_salt -
            3aq.NH43HSO42.M_salt - aq.NH4Cl.M_salt)/2, M_zero),
            max(SO4_extra - aq.NH4HSO4.M_salt - aq.Na2SO4.M_salt - aq.K2SO4.M_salt -
            aq.NaHSO4.M_salt - 2aq.NH43HSO42.M_salt, M_zero))
        aq.KCl.M_salt ~ min(max(K_extra - 2aq.K2SO4.M_salt, M_zero),
            max(Cl_extra - 2aq.CaCl2.M_salt - 2aq.MgCl2.M_salt - aq.NaCl.M_salt -
                aq.NH4Cl.M_salt, M_zero))
        aq.KHSO4.M_salt ~ min(max(K_extra - 2aq.K2SO4.M_salt - aq.KCl.M_salt, M_zero),
            max(SO4_extra - aq.K2SO4.M_salt  - aq.NH4HSO4.M_salt - aq.Na2SO4.M_salt -
            aq.NaHSO4.M_salt - 2aq.NH43HSO42.M_salt - aq.NH42SO4.M_salt, M_zero))
        aq.MgSO4.M_salt ~ min(max(Mg_extra - aq.MgCl2.M_salt - aq.MgNO32.M_salt, M_zero),
            max(SO4_extra - aq.KHSO4.M_salt- aq.K2SO4.M_salt  - aq.NH4HSO4.M_salt -
            aq.Na2SO4.M_salt - aq.NaHSO4.M_salt - 2aq.NH43HSO42.M_salt -
            aq.NH42SO4.M_salt, M_zero))
        aq.KNO3.M_salt ~ min(
            max(K_extra - 2aq.K2SO4.M_salt - aq.KCl.M_salt - aq.KHSO4.M_salt, M_zero),
            max(NO3_extra - aq.NaNO3.M_salt - 2aq.CaNO32.M_salt - 2aq.MgNO32.M_salt -
            aq.NH4NO3.M_salt, M_zero))
       aq.CaSO4.M_salt ~ min(max(Ca_extra - aq.CaCl2.M_salt - aq.CaNO32.M_salt, M_zero),
            max(SO4_extra - aq.MgSO4.M_salt - aq.KHSO4.M_salt- aq.K2SO4.M_salt  -
            aq.NH4HSO4.M_salt - aq.Na2SO4.M_salt -
            aq.NaHSO4.M_salt - 2aq.NH43HSO42.M_salt - aq.NH42SO4.M_salt, M_zero))

     #aq.CaCl2.M_salt ~ 0
#aq.MgCl2.M_salt ~ 0
   #     aq.NH4HSO4.M_salt ~ 0
   # aq.CaNO32.M_salt ~ 0
    #aq.Na2SO4.M_salt ~ 0
    #aq.NaHSO4.M_salt ~ 0
    #aq.MgNO32.M_salt ~ 0
    #aq.NH4NO3.M_salt ~ 0
    #aq.NH43HSO42.M_salt ~ 0
    #aq.NaNO3.M_salt ~ 0
    #aq.NaCl.M_salt ~ 0
    #aq.NH4Cl.M_salt ~ 0
   # aq.NH42SO4.M_salt ~ 0
   # aq.K2SO4.M_salt ~ 0
    #aq.KCl.M_salt ~ 0
#aq.KHSO4.M_salt ~ 0
      #  aq.MgSO4.M_salt ~ 0
       # aq.KNO3.M_salt ~ 0
    #    aq.CaSO4.M_salt ~ 0

        # These salts don't exist in solid form at atmospheric conditions.
        # We will use them to absorb any remaining extra mass.
        aq.HHSO4.M_salt ~ 0
        aq.H2SO4.M_salt ~ 0
        aq.HNO3.M_salt ~ 0
        aq.HCl.M_salt ~ 0
        aq.NH3_dissociated.M_salt ~ 0
        aq.H2O_dissociated.M_salt ~ 0
        aq.HSO4_dissociated.M_salt ~ 0

        # This is what's left over after allocating as much mass as possible to the salts.
        SO4_resid ~ SO4_extra - aq.MgSO4.M_salt - aq.KHSO4.M_salt-
            aq.K2SO4.M_salt  - aq.NH4HSO4.M_salt - aq.Na2SO4.M_salt - aq.NaHSO4.M_salt -
            2aq.NH43HSO42.M_salt - aq.NH42SO4.M_salt - aq.CaSO4.M_salt
        NO3_resid ~ NO3_extra - aq.NaNO3.M_salt - 2aq.CaNO32.M_salt -
            2aq.MgNO32.M_salt - aq.NH4NO3.M_salt -aq.KNO3.M_salt
         NH_resid ~NH_extra - aq.NH4HSO4.M_salt - aq.NH4NO3.M_salt -
            3aq.NH43HSO42.M_salt - aq.NH4Cl.M_salt - 2aq.NH42SO4.M_salt
        Cl_resid ~ Cl_extra - 2aq.CaCl2.M_salt - 2aq.MgCl2.M_salt - aq.NaCl.M_salt -
                aq.NH4Cl.M_salt - aq.KCl.M_salt
    end
end

@doc """
    Isorropia(kwargs...)

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

end
