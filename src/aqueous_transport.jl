export AqueousDiffusionReaction, MassTransportLimitation, DropletMassBalance

"""
    AqueousDiffusionReaction(; name=:AqueousDiffusionReaction)

Calculate the correction factor Q for aqueous-phase reaction rate when
diffusion within the droplet limits the reaction.

The factor Q accounts for concentration gradients within the droplet
due to fast reactions near the surface depleting reactants.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.110
"""
@component function AqueousDiffusionReaction(; name = :AqueousDiffusionReaction)
    @parameters begin
        R_p = 1.0e-5, [description = "Droplet radius", unit = u"m"]
        k_rxn = 1.0, [description = "First-order reaction rate constant", unit = u"s^-1"]
        D_aq = 1.0e-9, [description = "Aqueous-phase diffusivity", unit = u"m^2/s"]
    end

    @variables begin
        q(t),
        [
            description = "Dimensionless reaction-diffusion parameter (dimensionless)", unit = u"1"]
        Q(t),
        [description = "Aqueous diffusion correction factor (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 12.105: q = R_p √(k/D_aq)
        q ~ R_p * sqrt(k_rxn / D_aq),

        # Eq. 12.110: Q = 3(coth(q)/q - 1/q²)
        # Note: For numerical stability, this formula approaches 1 as q → 0
        # coth(q) = (exp(q) + exp(-q))/(exp(q) - exp(-q))
        Q ~ 3 * ((exp(q) + exp(-q)) / (exp(q) - exp(-q)) / q - 1 / q^2)
    ]

    return System(eqs, t; name)
end

"""
    MassTransportLimitation(; name=:MassTransportLimitation)

Evaluate mass transport limitation criteria for aqueous-phase chemistry.

Computes the limiting values of k₁H* below which mass transport does not
limit the aqueous-phase reaction rate (within 10%).

Reference: Seinfeld & Pandis (2006) Chapter 12, Eqs. 12.85, 12.86, 12.93
"""
@component function MassTransportLimitation(; name = :MassTransportLimitation)
    @constants begin
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
        ε = 0.1, [description = "Limitation threshold (10%) (dimensionless)", unit = u"1"]
    end

    @parameters begin
        R_p = 1.0e-5, [description = "Droplet radius", unit = u"m"]
        D_g = 2.0e-5, [description = "Gas-phase diffusivity", unit = u"m^2/s"]
        D_aq = 1.0e-9, [description = "Aqueous-phase diffusivity", unit = u"m^2/s"]
        α = 1.0, [description = "Accommodation coefficient (dimensionless)", unit = u"1"]
        M_A = 0.029, [description = "Molecular weight of species A", unit = u"kg/mol"]
        T = 298.15, [description = "Temperature", unit = u"K"]
    end

    @variables begin
        k1H_gas_limit(t),
        [description = "Gas-phase diffusion limitation threshold for k₁H*",
            unit = u"mol/(m^3*Pa*s)"]
        k1_aq_limit(t),
        [
            description = "Aqueous-phase diffusion limitation threshold for k₁", unit = u"s^-1"]
        k1H_interface_limit(t),
        [
            description = "Interfacial limitation threshold for k₁H*", unit = u"mol/(m^3*Pa*s)"]
    end

    eqs = [
        # Eq. 12.85: k_1 H*_A ≤ ε (3D_g)/(R_p² RT)
        k1H_gas_limit ~ ε * (3 * D_g) / (R_p^2 * R * T),

        # Eq. 12.86: k_1 ≤ ε (π²D_aq)/R_p²
        k1_aq_limit ~ ε * (π^2 * D_aq) / R_p^2,

        # Eq. 12.93: k_1 H*_A ≤ ε (3α)/(R_p √(2πM_A RT))
        k1H_interface_limit ~ ε * (3 * α) / (R_p * sqrt(2 * π * M_A * R * T))
    ]

    return System(eqs, t; name)
end

"""
    DropletMassBalance(; name=:DropletMassBalance)

Coupled gas-phase and aqueous-phase mass balance for cloud droplet chemistry.

This model couples:

  - Gas-phase partial pressure evolution
  - Aqueous-phase concentration evolution
  - Mass transfer between phases
  - Aqueous-phase chemical reaction

Reference: Seinfeld & Pandis (2006) Chapter 12, Eqs. 12.122-12.123
"""
@component function DropletMassBalance(; name = :DropletMassBalance)
    @constants begin
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
    end

    @parameters begin
        k_mt = 1.0, [description = "Mass transfer coefficient", unit = u"s^-1"]
        w_L = 1.0e-6,
        [description = "Liquid water volume fraction (dimensionless)", unit = u"1"]
        H_star = 1.0e5,
        [description = "Effective Henry's law coefficient", unit = u"mol/(m^3*Pa)"]
        T = 298.15, [description = "Temperature", unit = u"K"]
        Q = 1.0,
        [description = "Aqueous diffusion correction factor (dimensionless)", unit = u"1"]
        R_aq = 0.0, [description = "Aqueous-phase reaction rate", unit = u"mol/(m^3*s)"]
    end

    @variables begin
        p(t), [description = "Bulk gas-phase partial pressure", unit = u"Pa"]
        C_aq(t), [description = "Aqueous-phase concentration", unit = u"mol/m^3"]
    end

    eqs = [
        # Eq. 12.122: dp/dt = -k_mt w_L p + (1/H*) k_mt C_aq w_L
        D(p) ~ -k_mt * w_L * p + (1 / H_star) * k_mt * C_aq * w_L,

        # Eq. 12.123: dC_aq/dt = (k_mt/RT)p - (k_mt)/(H* RT) C_aq - Q R_aq
        D(C_aq) ~ (k_mt / (R * T)) * p - (k_mt / (H_star * R * T)) * C_aq - Q * R_aq
    ]

    return System(eqs, t; name)
end
