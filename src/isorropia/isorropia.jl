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

function salt_group(sys, group, var)
    salts = Dict(
        :NH4 => [:NH4NO3, :NH4Cl, :NH4HSO4, :NH42SO4, :NH43HSO42],
        :Na => [:NaCl, :Na2SO4, :NaNO3, :NaHSO4],
        :Ca => [:CaNO32, :CaCl2, :CaSO4],
        :K => [:KHSO4, :K2SO4, :KNO3, :KCl],
        :Mg => [:MgSO4, :MgNO32, :MgCl2],
        :Cl => [:NaCl, :KCl, :MgCl2, :CaCl2, :NH4Cl],
        :NO3 => [:NaNO3, :KNO3, :MgNO32, :CaNO32, :NH4NO3],
        :SO4 => [:Na2SO4, :K2SO4, :MgSO4, :CaSO4, :NH42SO4],
        :HSO4 => [:NH43HSO42, :KHSO4, :NaHSO4, :NH4HSO4, :NH43HSO42],
    )[group]
    return [reduce(getproperty, [sys, s, var]) for s in salts]
end

function salt_group_ν(group)
    return ν = Dict(
        :NH4 => :ν_cation,
        :Na => :ν_cation,
        :Ca => :ν_cation,
        :K => :ν_cation,
        :Mg => :ν_cation,
        :Cl => :ν_anion,
        :NO3 => :ν_anion,
        :SO4 => :ν_anion,
        :HSO4 => :ν_anion,
        :HSO42 => :ν_anion,
    )[group]
end

function resid_moles(sys, val, group, var, prior_vars, quantity, q_zero)
    ν = reduce(getproperty, [sys, var, salt_group_ν(group)])

    priors = [reduce(getproperty, [sys, v, quantity]) for v in prior_vars]
    prior_νs = [reduce(getproperty, [sys, v, salt_group_ν(group)]) for v in prior_vars]

    return max(val / ν - sum(priors .* prior_νs, init = 0), q_zero)
end

function min_resid(sys, var, cat, an, cat_val, an_val, cat_priors, an_priors, q_zero)
    salt = Dict(
        (:Ca, :NO3) => :CaNO32,
        (:Ca, :Cl) => :CaCl2,
        (:Ca, :SO4) => :CaSO4,
        (:K, :SO4) => :K2SO4,
        (:K, :NO3) => :KNO3,
        (:K, :Cl) => :KCl,
        (:Mg, :SO4) => :MgSO4,
        (:Mg, :NO3) => :MgNO32,
        (:Mg, :Cl) => :MgCl2,
        (:Na, :Cl) => :NaCl,
        (:Na, :SO4) => :Na2SO4,
        (:Na, :NO3) => :NaNO3,
        (:NH4, :HSO4) => :NH4HSO4,
        (:NH4, :NO3) => :NH4NO3,
        (:NH4, :Cl) => :NH4Cl,
        (:H, :HSO4) => :HHSO4,
        (:K, :HSO4) => :KHSO4,
        (:NH4, :SO4) => :NH42SO4,
        (:NH4, :HSO4) => :NH4HSO4,
        (:NH4, :HSO42) => :NH43HSO42,
        (:Na, :HSO4) => :NaHSO4,
    )[(cat, an)]
    return min(
        resid_moles(sys, cat_val, cat, salt, cat_priors, var, q_zero),
        resid_moles(sys, an_val, an, salt, an_priors, var, q_zero),
    )
end
function min_resid(
        sys,
        var,
        cat,
        an,
        cat_val,
        an_val,
        cat_priors,
        an_priors,
        W,
        m_one,
        q_zero,
    )
    M = min_resid(sys, var, cat, an, cat_val, an_val, cat_priors, an_priors, q_zero)
    return log(max(M / W / m_one, 1.0e-20))
end

@component function Species(; name = :Species, M_total = nothing)
    @variables begin
        total(t) = M_total, [unit = u"mol/m^3", description = "Total concentration"]
    end
    eqs = [
        D(total) ~ 0,
    ]
    return System(eqs, t; name)
end

@component function Isorropia(; name = :Isorropia)
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
        w_inf = 1.0e30, [unit = u"kg/m^3"]
    end

    aq = Aqueous(; name = :aq)
    g = Gases(; name = :g)
    eq = EquilibriumConstants(; name = :eq)

    NH = Species(; name = :NH, M_total = 2.0e-7)
    Na = Species(; name = :Na, M_total = 4.0e-10)
    Ca = Species(; name = :Ca, M_total = 1.0e-8)
    K = Species(; name = :K, M_total = 8.4e-9)
    Mg = Species(; name = :Mg, M_total = 4.1e-10)
    Cl = Species(; name = :Cl, M_total = 2.7e-10)
    NO3 = Species(; name = :NO3, M_total = 3.2e-8)
    SO4 = Species(; name = :SO4, M_total = 1.0e-7)

    @variables begin
        R_1(t), [description = "Total sulfate ratio from (Section 3.1) (dimensionless)", guess = 2.0]
        R_2(t), [description = "Crustal species and sodium ratio (Section 3.1) (dimensionless)", guess = 0.5]
        R_3(t), [description = "Crustal species ratio (Section 3.1) (dimensionless)", guess = 0.5]

        type1(t), [description = "Sulfate rich (free acid) aerosol type (Table 3) (dimensionless)", guess = 0.0]
        type2(t), [description = "Sulfate rich aerosol type (Table 3) (dimensionless)", guess = 0.0]
        type3(t), [description = "Sulfate poor, crustal & sodium poor aerosol (Table 3) (dimensionless)", guess = 1.0]
        type4(t), [description = "Sulfate poor, crustal & sodium rich, crustal poor (dimensionless)", guess = 0.0]
        type5(t), [description = "Sulfate poor, crustal & sodium rich, crustal rich (dimensionless)", guess = 0.0]

        NH_precip(t), [unit = u"mol/m^3", description = "Precipitated NH", guess = 0.0]
        Na_precip(t), [unit = u"mol/m^3", description = "Precipitated Na", guess = 0.0]
        Ca_precip(t), [unit = u"mol/m^3", description = "Precipitated Ca", guess = 0.0]
        K_precip(t), [unit = u"mol/m^3", description = "Precipitated K", guess = 0.0]
        Mg_precip(t), [unit = u"mol/m^3", description = "Precipitated Mg", guess = 0.0]
        Cl_precip(t), [unit = u"mol/m^3", description = "Precipitated Cl", guess = 0.0]
        NO3_precip(t), [unit = u"mol/m^3", description = "Precipitated NO3", guess = 0.0]
        SO4_precip(t), [unit = u"mol/m^3", description = "Precipitated SO4", guess = 0.0]
    end

    eqs = [
        # Connect subsystem variables to parent parameters
        aq.T ~ T,
        aq.RH ~ RH,
        g.T ~ T,
        eq.T ~ T,

        # Reactions based on information in Table 2 of Fountoukis and Nenes (2007).
        # The left-hand side of the reaction equation is the equilibrium constant and the
        # right-hand side is the ratio of the product and reactant activities.
        #eq.r1.logK_eq ~ aq.CaNO32.loga
        #eq.r2.logK_eq ~ aq.CaCl2.loga
        #eq.r3.logK_eq ~ aq.CaSO4.loga + 2log(RH)
        eq.r4.logK_eq ~ aq.K2SO4.loga,
        # eq.r5.logK_eq ~ aq.KHSO4.loga
        # eq.r6.logK_eq ~ aq.KNO3.loga
        # eq.r7.logK_eq ~ aq.KCl.loga
        # eq.r8.logK_eq ~ aq.MgSO4.loga
        # eq.r9.logK_eq ~ aq.MgNO32.loga
        # eq.r10.logK_eq ~ aq.MgCl2.loga
        eq.r11.logK_eq ~ aq.HSO4_dissociated.loga - aq.H2SO4.loga, # H2SO4 fully dissociates, so the concentration of HSO4 is the same as the concentration of HHSO4.
        eq.r12.logK_eq ~ aq.NH3.logm - g.NH3.logp,
        eq.r13.logK_eq ~ aq.NH3_dissociated.loga - aq.NH3.logm + log(RH),
        eq.r14.logK_eq ~ aq.HNO3.loga - g.HNO3.logp,
        eq.r15.logK_eq ~ aq.HNO3_aq.logm - g.HNO3.logp,
        eq.r16.logK_eq ~ aq.HCl.loga - g.HCl.logp,
        eq.r17.logK_eq ~ aq.HCl_aq.logm - g.HCl.logp,
        eq.r18.logK_eq ~ aq.H2O_dissociated.loga - log(RH),
        eq.r19.logK_eq ~ aq.Na2SO4.loga,
        # eq.r20.logK_eq ~ aq.NH42SO4.loga
        # eq.r22.logK_eq ~ aq.NaNO3.loga
        # eq.r23.logK_eq ~ aq.NaCl.loga
        # eq.r24.logK_eq ~ aq.NaHSO4.loga
        # eq.r26.logK_eq ~ aq.NH4HSO4.loga
        # eq.r27.logK_eq ~ aq.NH43HSO42.loga

        # These equations would result in an over-specified equilibrium.
        #eq.r21.logK_eq ~ log(g.NH3.p / p_one) + log(g.HCl.p / p_one)
        # eq.r25.logK_eq ~ log(g.NH3.p / p_one) + log(g.HNO3.p / p_one)

        ## Mass balance
        NH.total ~ g.NH3.M + aq.NH3.M + aq.NH3_dissociated.M, #+ NH.precip
        Cl.total ~ g.HCl.M + aq.HCl.M + aq.HCl_aq.M, #+ Cl.precip
        NO3.total ~ g.HNO3.M + aq.HNO3.M + aq.HNO3_aq.M, #+ NO3.precip
        SO4.total ~ aq.HSO4_dissociated.M + aq.H2SO4.M, #+ SO4.precip

        # D(*.total) ~ 0 equations are defined in Species

        # Aerosol types from Section 3.1 and Table 3
        R_1 ~ (NH.total + Ca.total + K.total + Mg.total + Na.total) / SO4.total,
        R_2 ~ (Ca.total + K.total + Mg.total + Na.total) / SO4.total,
        R_3 ~ (Ca.total + K.total + Mg.total) / SO4.total,
        type1 ~ 1 - (tanh((R_1 - 1) * 30) + 1) / 2,
        type2 ~ min((tanh((R_1 - 1) * 30) + 1) / 2, 1 - (tanh((R_1 - 2) * 30) + 1) / 2),
        type3 ~ min((tanh((R_1 - 2) * 30) + 1) / 2, 1 - (tanh((R_2 - 2) * 30) + 1) / 2),
        type4 ~ min(
            (tanh((R_1 - 2) * 30) + 1) / 2,
            (tanh((R_2 - 2) * 30) + 1) / 2,
            1 - (tanh((R_3 - 2) * 30) + 1) / 2,
        ),
        type5 ~ min(
            (tanh((R_1 - 2) * 30) + 1) / 2,
            (tanh((R_2 - 2) * 30) + 1) / 2,
            (tanh((R_3 - 2) * 30) + 1) / 2,
        ),

        # Aqueous mass balances
        # First fill in salts that are needed to guarantee space, or that are otherwise prioritized.
        aq.Na2SO4.M₀ ~
            min_resid(aq, :M, :Na, :SO4, Na.total, aq.HSO4_dissociated.M, [], [], M_zero), # Na preferentially forms Na2SO4 over other salts (section 3.2)
        aq.K2SO4.M₀ ~
            min_resid(aq, :M, :K, :SO4, K.total, aq.HSO4_dissociated.M, [], [:Na2SO4], M_zero), # K preferentially forms K2SO4 over other salts (section 3.2)
        aq.CaCl2.M ~ min_resid(aq, :M, :Ca, :Cl, Ca.total, aq.HCl.M, [], [], M_zero),
        aq.MgCl2.M ~ min_resid(aq, :M, :Mg, :Cl, Mg.total, aq.HCl.M, [], [:CaCl2], M_zero),
        aq.MgNO32.M ~
            min_resid(aq, :M, :Mg, :NO3, Mg.total, aq.HNO3.M, [:MgCl2], [], M_zero),
        aq.NH4HSO4.M ~
            min_resid(aq, :M, :NH4, :HSO4, aq.NH3_dissociated.M, aq.H2SO4.M, [], [], M_zero),
        aq.CaNO32.M ~
            min_resid(aq, :M, :Ca, :NO3, Ca.total, aq.HNO3.M, [:CaCl2], [:MgNO32], M_zero),
        aq.NaHSO4.M ~
            min_resid(aq, :M, :Na, :HSO4, Na.total, aq.H2SO4.M, [:Na2SO4], [:NH4HSO4], M_zero),
        # Next, fill in salts in order of increasing DRH unless otherwise noted.
        aq.NH4NO3.M ~ min_resid(
            aq,
            :M,
            :NH4,
            :NO3,
            aq.NH3_dissociated.M,
            aq.HNO3.M,
            [:NH4HSO4],
            [:MgNO32, :CaNO32],
            M_zero,
        ),
        aq.NH43HSO42.M ~ min_resid(
            aq,
            :M,
            :NH4,
            :HSO42,
            aq.NH3_dissociated.M,
            aq.H2SO4.M,
            [:NH4HSO4, :NH4NO3],
            [:NH4HSO4, :NaHSO4],
            M_zero,
        ),
        aq.NaNO3.M ~ min_resid(
            aq,
            :M,
            :Na,
            :NO3,
            Na.total,
            aq.HNO3.M,
            [:Na2SO4, :NaHSO4],
            [:MgNO32, :CaNO32, :NH4NO3],
            M_zero,
        ),
        aq.NaCl.M ~ min_resid(
            aq,
            :M,
            :Na,
            :Cl,
            Na.total,
            aq.HCl.M,
            [:Na2SO4, :NaHSO4, :NaNO3],
            [:CaCl2, :MgCl2],
            M_zero,
        ),
        aq.NH4Cl.M ~ min_resid(
            aq,
            :M,
            :NH4,
            :Cl,
            aq.NH3_dissociated.M,
            aq.HCl.M,
            [:NH4HSO4, :NH4NO3, :NH43HSO42],
            [:CaCl2, :MgCl2, :NaCl],
            M_zero,
        ),
        aq.NH42SO4.M ~ min_resid(
            aq,
            :M,
            :NH4,
            :SO4,
            aq.NH3_dissociated.M,
            aq.HSO4_dissociated.M,
            [:NH4HSO4, :NH4NO3, :NH43HSO42, :NH4Cl],
            [:Na2SO4, :K2SO4],
            M_zero,
        ),
        aq.KCl.M ~ min_resid(
            aq,
            :M,
            :K,
            :Cl,
            K.total,
            aq.HCl.M,
            [:K2SO4],
            [:CaCl2, :MgCl2, :NaCl, :NH4Cl],
            M_zero,
        ),
        aq.KHSO4.M ~ min_resid(
            aq,
            :M,
            :K,
            :HSO4,
            K.total,
            aq.H2SO4.M,
            [:K2SO4, :KCl],
            [:NH4HSO4, :NaHSO4, :NH43HSO42],
            M_zero,
        ),
        aq.MgSO4.M ~ min_resid(
            aq,
            :M,
            :Mg,
            :SO4,
            Mg.total,
            aq.HSO4_dissociated.M,
            [:MgCl2, :MgNO32],
            [:Na2SO4, :K2SO4, :NH42SO4],
            M_zero,
        ),
        aq.KNO3.M ~ min_resid(
            aq,
            :M,
            :K,
            :NO3,
            K.total,
            aq.HNO3.M,
            [:K2SO4, :KCl, :KHSO4],
            [:MgNO32, :CaNO32, :NH4NO3, :NaNO3],
            M_zero,
        ),
        aq.CaSO4.M ~ min_resid(
            aq,
            :M,
            :Ca,
            :SO4,
            Ca.total,
            aq.HSO4_dissociated.M,
            [:CaCl2, :CaNO32],
            [:Na2SO4, :K2SO4, :NH42SO4],
            M_zero,
        ),

        NH_precip ~
            aq.NH3_dissociated.M -
            sum(salt_group(aq, :NH4, :M) .* salt_group(aq, :NH4, salt_group_ν(:NH4))),
        Na_precip ~
            Na.total - sum(salt_group(aq, :Na, :M) .* salt_group(aq, :Na, salt_group_ν(:Na))),
        Ca_precip ~
            Ca.total - sum(salt_group(aq, :Ca, :M) .* salt_group(aq, :Ca, salt_group_ν(:Ca))),
        K_precip ~
            K.total - sum(salt_group(aq, :K, :M) .* salt_group(aq, :K, salt_group_ν(:K))),
        Mg_precip ~
            Mg.total - sum(salt_group(aq, :Mg, :M) .* salt_group(aq, :Mg, salt_group_ν(:Mg))),
        Cl_precip ~
            Cl.total - sum(salt_group(aq, :Cl, :M) .* salt_group(aq, :Cl, salt_group_ν(:Cl))),
        NO3_precip ~
            aq.HNO3.M -
            sum(salt_group(aq, :NO3, :M) .* salt_group(aq, :NO3, salt_group_ν(:NO3))),
        SO4_precip ~
            aq.HSO4_dissociated.M -
            sum(salt_group(aq, :SO4, :M) .* salt_group(aq, :SO4, salt_group_ν(:SO4))) +
            aq.H2SO4.M -
            sum(salt_group(aq, :HSO4, :M) .* salt_group(aq, :HSO4, salt_group_ν(:HSO4))),
    ]

    return System(eqs, t; systems = [aq, g, eq, NH, Na, Ca, K, Mg, Cl, NO3, SO4], name)
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
    :Na2SO4_s => 142.0421,
)

end
