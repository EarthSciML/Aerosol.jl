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

@mtkmodel Species begin
    @structural_parameters begin
        total_guess=nothing
    end
    @variables begin
        total(t), [unit = u"mol/m^3", description="Total concentration", guess=total_guess]
        eq(t), [unit = u"mol/kg", description="Concentration at equilibrium", guess=1]
        #! format: off
        extra(t), [unit = u"mol/m^3", description = "Extra concentration above equilibrium for allocation to salts"]
        #! format: on
        resid(t), [unit = u"mol/m^3", description = "Residual conc. after salt formation"]
        # aq(t), [unit = u"mol/m^3", description="Aqueous concentration"]
        # precip(t), [unit = u"mol/m^3", description="Precipitated (solid) concentration"]
        # aer(t), [unit=u"mol/m^3", description="Aerosol concentration"]
        # f_aer(t), [description = "Mole fraction aerosol (vs. gas)"]
    end
    # @equations begin
    #     aer ~ aq + precip
    #     f_aer ~ aer / total
    # end
end

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

        NH = Species(total_guess=2e-7)
        Na = Species(total_guess=4e-10)
        Ca = Species(total_guess=1e-8)
        K = Species(total_guess=8.4e-9)
        Mg = Species(total_guess=4.1e-10)
        Cl = Species(total_guess=2.7e-10)
        NO3 = Species(total_guess=3.2e-8)
        SO4 = Species(total_guess=1e-7)
    end
    @variables begin
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
        eq.r14.logK_eq ~ aq.HNO3.loga_aq - log(g.HNO3.p / p_one) - 16 # FIXME(CT): Added -16 to get reasonable results; not sure why this is needed.
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
        #eq.r11.logK_eq ~ 2log(aq.HSO4_dissociated.a / m_one) - aq.H2SO4.loga # H2SO4 fully dissociates, so the concentration of HSO4 is the same as the concentration of HHSO4.
        #eq.r21.logK_eq ~ log(g.NH3.p / p_one) + log(g.HCl.p / p_one)
        # eq.r25.logK_eq ~ log(g.NH3.p / p_one) + log(g.HNO3.p / p_one)

        # Aqueous equilibrium mass Balance
        NH.eq ~ aq.NH3_dissociated.m_eq + aq.NH3.m_eq
        Na.eq ~ aq.Na.m_eq
        Ca.eq ~ aq.Ca.m_eq
        K.eq ~ aq.K.m_eq
        Mg.eq ~ aq.Mg.m_eq
        Cl.eq ~ aq.HCl.m_eq + aq.HCl_aq.m_eq
        NO3.eq ~ aq.HNO3.m_eq + aq.HNO3_aq.m_eq
        SO4.eq ~ aq.HSO4_dissociated.m_eq + aq.HHSO4.m_eq

        NH.extra ~ NH.total - g.NH3.M - NH.eq * aq.W_eq
        Na.extra ~ Na.total - Na.eq * aq.W_eq
        Ca.extra ~ Ca.total - Ca.eq * aq.W_eq
        K.extra ~ K.total - K.eq * aq.W_eq
        Mg.extra ~ Mg.total - Mg.eq * aq.W_eq
        Cl.extra ~ Cl.total - g.HCl.M - Cl.eq * aq.W_eq
        NO3.extra ~ NO3.total - g.HNO3.M - NO3.eq * aq.W_eq
        SO4.extra ~ SO4.total - SO4.eq * aq.W_eq
        aq.W_eq ~ min(
            (NH.total - g.NH3.M) / NH.eq,
            Na.total / Na.eq,
            Ca.total / Ca.eq,
            K.total / K.eq,
            Mg.total / Mg.eq,
            (Cl.total - g.HCl.M) / Cl.eq,
            (NO3.total - g.HNO3.M) / NO3.eq,
            SO4.total / SO4.eq
        )

        D(NH.total) ~ 0.0
        D(Na.total) ~ 0.0
        D(Ca.total) ~ 0.0
        D(K.total) ~ 0.0
        D(Mg.total) ~ 0.0
        D(Cl.total) ~ 0.0
        D(NO3.total) ~ 0.0
        D(SO4.total) ~ 0.0

        # Aerosol types from Section 3.1 and Table 3
        R_1 ~ (NH.total + Ca.total + K.total + Mg.total + Na.total) / SO4.total
        R_2 ~ (Ca.total + K.total + Mg.total + Na.total) / SO4.total
        R_3 ~ (Ca.total + K.total + Mg.total) / SO4.total
        type1 ~ 1 - (tanh((R_1 - 1) * 30) + 1) / 2
        type2 ~ min((tanh((R_1 - 1) * 30) + 1) / 2, 1 - (tanh((R_1 - 2) * 30) + 1) / 2)
        type3 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, 1 - (tanh((R_2 - 2) * 30) + 1) / 2)
        type4 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, (tanh((R_2 - 2) * 30) + 1) / 2,
            1 - (tanh((R_3 - 2) * 30) + 1) / 2)
        type5 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, (tanh((R_2 - 2) * 30) + 1) / 2,
            1 - (tanh((R_3 - 2) * 30) + 1) / 2)

        # Solid mass balances
        # First fill in salts that are needed to guarantee space, or that are otherwise prioritized.
        aq.Na2SO4.M_salt ~ min(Na.extra/2, SO4.extra) # Na preferentially forms Na2SO4 over other salts (section 3.2)
        aq.K2SO4.M_salt ~ min(K.extra/2, max(SO4.extra - aq.Na2SO4.M_salt, M_zero)) # K preferentially forms K2SO4 over other salts (section 3.2)
        aq.CaCl2.M_salt ~ min(Ca.extra, Cl.extra/2)
         aq.MgCl2.M_salt ~ min(Mg.extra, Cl.extra/2 - aq.CaCl2.M_salt)
         aq.MgNO32.M_salt ~ min(Mg.extra - aq.MgCl2.M_salt, NO3.extra/2)
        aq.NH4HSO4.M_salt ~ min(NH.extra, max(SO4.extra - aq.Na2SO4.M_salt -
            aq.K2SO4.M_salt, M_zero))
         aq.CaNO32.M_salt ~ min(Ca.extra, max(NO3.extra/2 - aq.CaNO32.M_salt, M_zero))
        aq.NaHSO4.M_salt ~ min(max(Na.extra - 2aq.Na2SO4.M_salt, M_zero),
            max(SO4.extra - aq.NH4HSO4.M_salt - aq.Na2SO4.M_salt - aq.K2SO4.M_salt, M_zero))
        # Next, fill in salts in order of increasing DRH unless otherwise noted.
        aq.NH4NO3.M_salt ~ min(max(NH.extra - aq.NH4HSO4.M_salt, M_zero),
            max(NO3.extra - 2aq.CaNO32.M_salt - 2aq.MgNO32.M_salt, M_zero))
        aq.NH43HSO42.M_salt ~ min(
            max((NH.extra  - aq.NH4HSO4.M_salt- aq.NH4NO3.M_salt )/3, M_zero),
            max((SO4.extra - aq.Na2SO4.M_salt- aq.K2SO4.M_salt - aq.NH4HSO4.M_salt -
                aq.NH4HSO4.M_salt)/2, M_zero))
        aq.NaNO3.M_salt ~ min(max(Na.extra - 2aq.Na2SO4.M_salt - aq.NaHSO4.M_salt, M_zero),
            max(NO3.extra - 2aq.CaNO32.M_salt - 2aq.MgNO32.M_salt - aq.NH4NO3.M_salt, M_zero))
        aq.NaCl.M_salt ~ min(max(Na.extra - 2aq.Na2SO4.M_salt - aq.NaHSO4.M_salt -
            aq.NaNO3.M_salt, M_zero),
            max(Cl.extra - 2aq.CaCl2.M_salt - 2aq.MgCl2.M_salt, M_zero))
        aq.NH4Cl.M_salt ~ min(NH.extra - aq.NH4HSO4.M_salt - aq.NH4NO3.M_salt -
            3aq.NH43HSO42.M_salt,
            Cl.extra - 2aq.CaCl2.M_salt - 2aq.MgCl2.M_salt - aq.NaCl.M_salt)
        aq.NH42SO4.M_salt ~ min(max((NH.extra - aq.NH4HSO4.M_salt - aq.NH4NO3.M_salt -
            3aq.NH43HSO42.M_salt - aq.NH4Cl.M_salt)/2, M_zero),
            max(SO4.extra - aq.NH4HSO4.M_salt - aq.Na2SO4.M_salt - aq.K2SO4.M_salt -
            aq.NaHSO4.M_salt - 2aq.NH43HSO42.M_salt, M_zero))
        aq.KCl.M_salt ~ min(max(K.extra - 2aq.K2SO4.M_salt, M_zero),
            max(Cl.extra - 2aq.CaCl2.M_salt - 2aq.MgCl2.M_salt - aq.NaCl.M_salt -
                aq.NH4Cl.M_salt, M_zero))
        aq.KHSO4.M_salt ~ min(max(K.extra - 2aq.K2SO4.M_salt - aq.KCl.M_salt, M_zero),
            max(SO4.extra - aq.K2SO4.M_salt  - aq.NH4HSO4.M_salt - aq.Na2SO4.M_salt -
            aq.NaHSO4.M_salt - 2aq.NH43HSO42.M_salt - aq.NH42SO4.M_salt, M_zero))
        aq.MgSO4.M_salt ~ min(max(Mg.extra - aq.MgCl2.M_salt - aq.MgNO32.M_salt, M_zero),
            max(SO4.extra - aq.KHSO4.M_salt- aq.K2SO4.M_salt  - aq.NH4HSO4.M_salt -
            aq.Na2SO4.M_salt - aq.NaHSO4.M_salt - 2aq.NH43HSO42.M_salt -
            aq.NH42SO4.M_salt, M_zero))
        aq.KNO3.M_salt ~ min(
            max(K.extra - 2aq.K2SO4.M_salt - aq.KCl.M_salt - aq.KHSO4.M_salt, M_zero),
            max(NO3.extra - aq.NaNO3.M_salt - 2aq.CaNO32.M_salt - 2aq.MgNO32.M_salt -
            aq.NH4NO3.M_salt, M_zero))
       aq.CaSO4.M_salt ~ min(max(Ca.extra - aq.CaCl2.M_salt - aq.CaNO32.M_salt, M_zero),
            max(SO4.extra - aq.MgSO4.M_salt - aq.KHSO4.M_salt- aq.K2SO4.M_salt  -
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
        aq.HHSO4.M_salt ~ 0
        aq.H2SO4.M_salt ~ 0
        aq.HNO3.M_salt ~ 0
        aq.HCl.M_salt ~ 0
        aq.NH3_dissociated.M_salt ~ 0
        aq.H2O_dissociated.M_salt ~ 0
        aq.HSO4_dissociated.M_salt ~ 0

        # This is what's left over after allocating as much mass as possible to the salts.
        SO4.resid ~ SO4.extra - aq.MgSO4.M_salt - aq.KHSO4.M_salt-
            aq.K2SO4.M_salt  - aq.NH4HSO4.M_salt - aq.Na2SO4.M_salt - aq.NaHSO4.M_salt -
            2aq.NH43HSO42.M_salt - aq.NH42SO4.M_salt - aq.CaSO4.M_salt
        NO3.resid ~ NO3.extra - aq.NaNO3.M_salt - 2aq.CaNO32.M_salt -
            2aq.MgNO32.M_salt - aq.NH4NO3.M_salt -aq.KNO3.M_salt
         NH.resid ~NH.extra - aq.NH4HSO4.M_salt - aq.NH4NO3.M_salt -
            3aq.NH43HSO42.M_salt - aq.NH4Cl.M_salt - 2aq.NH42SO4.M_salt
        Cl.resid ~ Cl.extra - 2aq.CaCl2.M_salt - 2aq.MgCl2.M_salt - aq.NaCl.M_salt -
                aq.NH4Cl.M_salt - aq.KCl.M_salt
        Na.resid ~ 0
        Ca.resid ~ 0
        K.resid ~ 0
        Mg.resid ~ 0

        # # Calculate dissolved and solid totals.
        # NH.aq ~ sum([x.ν_cation * x.M_aq for x in [
        #         aq.NH4NO3, aq.NH4Cl, aq.NH4HSO4, aq.NH42SO4, aq.NH43HSO42]])
        # NH.precip ~ sum([x.ν_cation * x.M_precip for x in [
        #         aq.NH4NO3, aq.NH4Cl, aq.NH4HSO4, aq.NH42SO4, aq.NH43HSO42]])
        # Na.aq ~sum([x.ν_cation * x.M_aq for x in [
        #     aq.NaCl, aq.Na2SO4, aq.NaNO3, aq.NaHSO4]])
        # Na.precip ~sum([x.ν_cation * x.M_precip for x in [
        #     aq.NaCl, aq.Na2SO4, aq.NaNO3, aq.NaHSO4]])
        # Ca.aq ~sum([x.ν_cation * x.M_aq for x in [
        #     aq.CaNO32, aq.CaCl2, aq.CaSO4]])
        # Ca.precip ~sum([x.ν_cation * x.M_precip for x in [
        #     aq.CaNO32, aq.CaCl2, aq.CaSO4]])
        # K.aq ~sum([x.ν_cation * x.M_aq for x in [
        #     aq.KHSO4, aq.K2SO4, aq.KNO3, aq.KCl]])
        # K.precip ~sum([x.ν_cation * x.M_precip for x in [
        #     aq.KHSO4, aq.K2SO4, aq.KNO3, aq.KCl]])
        # Mg.aq ~sum([x.ν_cation * x.M_aq for x in [
        #     aq.MgSO4, aq.MgNO32, aq.MgCl2]])
        # Mg.precip~sum([x.ν_cation * x.M_precip for x in [
        #     aq.MgSO4, aq.MgNO32, aq.MgCl2]])
        # Cl.aq ~sum([x.ν_anion * x.M_aq for x in [
        #     aq.NaCl, aq.KCl, aq.MgCl2, aq.CaCl2, aq.NH4Cl]])
        # Cl.precip ~sum([x.ν_anion * x.M_precip for x in [
        #     aq.NaCl, aq.KCl, aq.MgCl2, aq.CaCl2, aq.NH4Cl]])
        # NO3.aq ~sum([x.ν_anion * x.M_aq for x in [
        #     aq.NaNO3, aq.KNO3, aq.MgNO32, aq.CaNO32, aq.NH4NO3]])
        # NO3.precip ~sum([x.ν_anion * x.M_precip for x in [
        #     aq.NaNO3, aq.KNO3, aq.MgNO32, aq.CaNO32, aq.NH4NO3]])
        # SO4.aq ~sum([x.ν_anion * x.M_aq for x in [
        #     aq.Na2SO4, aq.K2SO4, aq.MgSO4, aq.CaSO4,aq.NH42SO4, aq.NH43HSO42,
        # aq.KHSO4, aq.NaHSO4, aq.NH4HSO4, aq.NH43HSO42]])
        # SO4.precip ~sum([x.ν_anion * x.M_precip for x in [
        #     aq.Na2SO4, aq.K2SO4, aq.MgSO4, aq.CaSO4,aq.NH42SO4, aq.NH43HSO42,
        # aq.KHSO4, aq.NaHSO4, aq.NH4HSO4, aq.NH43HSO42]])
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
