export GasDiffusionTimescale, AqueousDiffusionTimescale, InterfacialTimescale
export ReactionTimescale, SolidEquilibrationTimescale, AqueousEquilibrationTimescale

"""
    GasDiffusionTimescale(; name=:GasDiffusionTimescale)

Characteristic timescale for gas-phase diffusion to a particle.

This is the time required for a gas-phase concentration perturbation
to diffuse to the particle surface.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.49
"""
@component function GasDiffusionTimescale(; name = :GasDiffusionTimescale)
    @parameters begin
        R_p = 1.0e-5, [description = "Particle radius", unit = u"m"]
        D_g = 2.0e-5, [description = "Gas-phase diffusivity", unit = u"m^2/s"]
    end

    @variables begin
        τ_dg(t), [description = "Gas-phase diffusion timescale", unit = u"s"]
    end

    eqs = [
    # Eq. 12.49: τ_dg = R_p²/(4D_g)
        τ_dg ~ R_p^2 / (4 * D_g),
    ]

    return System(eqs, t; name)
end

"""
    AqueousDiffusionTimescale(; name=:AqueousDiffusionTimescale)

Characteristic timescale for aqueous-phase diffusion within a droplet.

This is the time for a dissolved species to diffuse throughout the droplet
interior.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.75
"""
@component function AqueousDiffusionTimescale(; name = :AqueousDiffusionTimescale)
    @parameters begin
        R_p = 1.0e-5, [description = "Droplet radius", unit = u"m"]
        D_aq = 1.0e-9, [description = "Aqueous-phase diffusivity", unit = u"m^2/s"]
    end

    @variables begin
        τ_da(t), [description = "Aqueous-phase diffusion timescale", unit = u"s"]
    end

    eqs = [
    # Eq. 12.75: τ_da = R_p²/(π²D_aq)
        τ_da ~ R_p^2 / (π^2 * D_aq),
    ]

    return System(eqs, t; name)
end

"""
    InterfacialTimescale(; name=:InterfacialTimescale)

Characteristic timescale to achieve equilibrium at the gas-particle interface.

For very soluble gases, the timescale is controlled by interfacial transport.
For low solubility gases, the timescale is controlled by aqueous diffusion.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eqs. 12.61-12.62
"""
@component function InterfacialTimescale(; name = :InterfacialTimescale)
    @constants begin
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
    end

    @parameters begin
        R_p = 1.0e-5, [description = "Droplet radius", unit = u"m"]
        α = 1.0, [description = "Accommodation coefficient (dimensionless)", unit = u"1"]
        M_A = 0.029, [description = "Molecular weight of species A", unit = u"kg/mol"]
        T = 298.15, [description = "Temperature", unit = u"K"]
        H_star = 1.0e5,
        [description = "Effective Henry's law coefficient", unit = u"mol/(m^3*Pa)"]
        D_aq = 1.0e-9, [description = "Aqueous-phase diffusivity", unit = u"m^2/s"]
    end

    @variables begin
        τ_p_soluble(t),
        [description = "Interfacial timescale for very soluble gases", unit = u"s"]
        τ_p_insoluble(t),
        [description = "Interfacial timescale for low solubility gases", unit = u"s"]
    end

    eqs = [
        # Eq. 12.61: τ_p ≈ (R_p H*_A √(2πM_A RT))/(3α) for very soluble gases
        τ_p_soluble ~ (R_p * H_star * sqrt(2 * π * M_A * R * T)) / (3 * α),

        # Eq. 12.62: τ_p ≈ R_p²/(π²D_aq) for low solubility gases
        τ_p_insoluble ~ R_p^2 / (π^2 * D_aq)
    ]

    return System(eqs, t; name)
end

"""
    ReactionTimescale(; name=:ReactionTimescale)

Characteristic timescale for aqueous-phase chemical reactions.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eqs. 12.76-12.78
"""
@component function ReactionTimescale(; name = :ReactionTimescale)
    @constants begin
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
    end

    @parameters begin
        k_rxn = 1.0, [description = "First-order reaction rate constant", unit = u"s^-1"]
        H_star = 1.0e5,
        [description = "Effective Henry's law coefficient", unit = u"mol/(m^3*Pa)"]
        T = 298.15, [description = "Temperature", unit = u"K"]
    end

    @variables begin
        τ_ra(t), [description = "Aqueous-phase reaction timescale", unit = u"s"]
        τ_rg(t), [description = "Gas-phase reference reaction timescale", unit = u"s"]
    end

    eqs = [
        # Eq. 12.76: τ_ra = 1/k (for first-order reaction)
        τ_ra ~ 1 / k_rxn,

        # Eq. 12.78: τ_rg = τ_ra/(H*RT)
        τ_rg ~ τ_ra / (H_star * R * T)
    ]

    return System(eqs, t; name)
end

"""
    SolidEquilibrationTimescale(; name=:SolidEquilibrationTimescale)

Characteristic timescale for a gas to equilibrate with solid aerosol particles.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eqs. 12.135, 12.139
"""
@component function SolidEquilibrationTimescale(; name = :SolidEquilibrationTimescale)
    @parameters begin
        R_p = 1.0e-7, [description = "Particle radius", unit = u"m"]
        ρ_p = 1000.0, [description = "Particle density", unit = u"kg/m^3"]
        D_A = 2.0e-5, [description = "Gas-phase diffusivity", unit = u"m^2/s"]
        m_p = 1.0e-8, [description = "Aerosol mass concentration", unit = u"kg/m^3"]
        f_Kn = 1.0,
        [
            description = "Transition regime correction factor f(Kn,α) (dimensionless)", unit = u"1"]
        N = 1.0e9, [description = "Particle number concentration", unit = u"m^-3"]
    end

    @variables begin
        τ_s(t), [description = "Solid particle equilibration timescale", unit = u"s"]
        τ_s_alt(t), [description = "Alternative form using N and R_p", unit = u"s"]
    end

    eqs = [
        # Eq. 12.135: τ_s = (ρ_p R_p²)/(3D_A m_p f(Kn,α))
        τ_s ~ (ρ_p * R_p^2) / (3 * D_A * m_p * f_Kn),

        # Eq. 12.139: τ_s ≈ 1/(4πN R̄_p D_A f̄)
        τ_s_alt ~ 1 / (4 * π * N * R_p * D_A * f_Kn)
    ]

    return System(eqs, t; name)
end

"""
    AqueousEquilibrationTimescale(; name=:AqueousEquilibrationTimescale)

Characteristic timescale for a gas to equilibrate with aqueous aerosol particles.

This accounts for the additional buffering effect of Henry's law dissolution.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.147
"""
@component function AqueousEquilibrationTimescale(; name = :AqueousEquilibrationTimescale)
    @parameters begin
        m_w = 1.0e-8, [description = "Aerosol water mass concentration", unit = u"kg/m^3"]
        K_A = 1.0e-3, [description = "Equilibrium constant", unit = u"kg/m^3"]
        τ_s = 1.0, [description = "Solid equilibration timescale", unit = u"s"]
    end

    @variables begin
        τ_a(t), [description = "Aqueous equilibration timescale", unit = u"s"]
    end

    eqs = [
    # Eq. 12.147: τ_a = (m_w/K_A) τ_s
        τ_a ~ (m_w / K_A) * τ_s,
    ]

    return System(eqs, t; name)
end
