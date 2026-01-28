export MassTransfer, MeanMolecularSpeed, KnudsenNumber, MeanFreePath
export FuchsSutugin, Dahneke, MaxwellianFlux, ParticleGrowthRate
export MassTransferCoefficient, UptakeCoefficient

"""
    MeanMolecularSpeed(; name=:MeanMolecularSpeed)

Calculate the mean molecular speed of gas molecules using kinetic theory.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.24
"""
@component function MeanMolecularSpeed(; name=:MeanMolecularSpeed)
    @constants begin
        k_B = 1.380649e-23, [description = "Boltzmann constant", unit = u"J/K"]
    end

    @parameters begin
        T = 298.15, [description = "Temperature", unit = u"K"]
        M_A = 0.029, [description = "Molecular weight of species A", unit = u"kg/mol"]
    end

    @constants begin
        N_Av = 6.02214076e23, [description = "Avogadro's number", unit = u"mol^-1"]
    end

    @variables begin
        c_bar(t), [description = "Mean molecular speed", unit = u"m/s"]
    end

    eqs = [
        # Eq. 12.24: c̄_A = √(8kT/(πm_A)) where m_A = M_A/N_Av
        c_bar ~ sqrt(8 * k_B * T * N_Av / (π * M_A)),
    ]

    return System(eqs, t; name)
end

"""
    MeanFreePath(; name=:MeanFreePath)

Calculate the mean free path of gas molecules in air.

Reference: Seinfeld & Pandis (2006) Chapter 12, Section 12.1.3
"""
@component function MeanFreePath(; name=:MeanFreePath)
    @parameters begin
        T = 298.15, [description = "Temperature", unit = u"K"]
        P = 101325.0, [description = "Pressure", unit = u"Pa"]
        M_A = 0.029, [description = "Molecular weight of species A", unit = u"kg/mol"]
        M_air = 0.029, [description = "Molecular weight of air", unit = u"kg/mol"]
        μ = 1.81e-5, [description = "Dynamic viscosity of air", unit = u"Pa*s"]
    end

    @constants begin
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
    end

    @variables begin
        λ(t), [description = "Mean free path", unit = u"m"]
        ρ_air(t), [description = "Air density", unit = u"kg/m^3"]
    end

    eqs = [
        # Air density from ideal gas law
        ρ_air ~ P * M_air / (R * T),
        # Mean free path: λ = 2μ/(ρ_air * c̄) where c̄ = √(8RT/(πM))
        # Simplified form commonly used
        λ ~ 2 * μ / (ρ_air * sqrt(8 * R * T / (π * M_air))),
    ]

    return System(eqs, t; name)
end

"""
    KnudsenNumber(; name=:KnudsenNumber)

Calculate the Knudsen number, the ratio of mean free path to particle radius.

The Knudsen number determines the transport regime:
- Kn ≪ 1: Continuum regime (Fick's law applies)
- Kn ≫ 1: Kinetic regime (molecular kinetics applies)
- Kn ~ 1: Transition regime (interpolation formulas needed)

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.31
"""
@component function KnudsenNumber(; name=:KnudsenNumber)
    @parameters begin
        λ = 6.5e-8, [description = "Mean free path", unit = u"m"]
        R_p = 1.0e-7, [description = "Particle radius", unit = u"m"]
    end

    @variables begin
        Kn(t), [description = "Knudsen number (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 12.31: Kn = λ/R_p
        Kn ~ λ / R_p,
    ]

    return System(eqs, t; name)
end

"""
    FuchsSutugin(; name=:FuchsSutugin)

Fuchs-Sutugin transition regime correction for mass transfer.

This is the most widely used formula for transition regime mass transfer due to its
simplicity and accuracy across the full Knudsen number range.

Returns the ratio J/J_c where J_c is the continuum regime flux.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.43 (with accommodation coefficient)
"""
@component function FuchsSutugin(; name=:FuchsSutugin)
    @parameters begin
        Kn = 1.0, [description = "Knudsen number (dimensionless)", unit = u"1"]
        α = 1.0, [description = "Accommodation coefficient (dimensionless)", unit = u"1"]
    end

    @variables begin
        f_FS(t), [description = "Fuchs-Sutugin correction factor J/J_c (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 12.43: J/J_c = 0.75α(1 + Kn) / (Kn² + Kn + 0.283Kn*α + 0.75α)
        f_FS ~ 0.75 * α * (1 + Kn) / (Kn^2 + Kn + 0.283 * Kn * α + 0.75 * α),
    ]

    return System(eqs, t; name)
end

"""
    Dahneke(; name=:Dahneke)

Dahneke transition regime correction for mass transfer.

An alternative interpolation formula for the transition regime.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.42
"""
@component function Dahneke(; name=:Dahneke)
    @parameters begin
        Kn = 1.0, [description = "Knudsen number (dimensionless)", unit = u"1"]
        α = 1.0, [description = "Accommodation coefficient (dimensionless)", unit = u"1"]
    end

    @variables begin
        f_D(t), [description = "Dahneke correction factor J/J_c (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 12.42: J/J_c = (1 + Kn) / (1 + 2Kn(1 + Kn)/α)
        f_D ~ (1 + Kn) / (1 + 2 * Kn * (1 + Kn) / α),
    ]

    return System(eqs, t; name)
end

"""
    MaxwellianFlux(; name=:MaxwellianFlux)

Calculate the Maxwellian (continuum regime) molar flux to a particle.

This is the steady-state diffusive flux in the continuum regime (Kn ≪ 1).

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.12
"""
@component function MaxwellianFlux(; name=:MaxwellianFlux)
    @parameters begin
        R_p = 1.0e-7, [description = "Particle radius", unit = u"m"]
        D_g = 2.0e-5, [description = "Gas-phase diffusivity", unit = u"m^2/s"]
        c_inf = 1.0e-6, [description = "Bulk gas-phase concentration", unit = u"mol/m^3"]
        c_s = 0.0, [description = "Surface concentration", unit = u"mol/m^3"]
    end

    @variables begin
        J_c(t), [description = "Continuum regime molar flow rate", unit = u"mol/s"]
    end

    eqs = [
        # Eq. 12.12: J_c = 4πR_p D_g (c_∞ - c_s)
        J_c ~ 4 * π * R_p * D_g * (c_inf - c_s),
    ]

    return System(eqs, t; name)
end

"""
    ParticleGrowthRate(; name=:ParticleGrowthRate)

Calculate the rate of particle radius change due to condensation/evaporation.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.14
"""
@component function ParticleGrowthRate(; name=:ParticleGrowthRate)
    @parameters begin
        D_g = 2.0e-5, [description = "Gas-phase diffusivity", unit = u"m^2/s"]
        M_A = 0.029, [description = "Molecular weight of species A", unit = u"kg/mol"]
        ρ_p = 1000.0, [description = "Particle density", unit = u"kg/m^3"]
        c_inf = 1.0e-6, [description = "Bulk gas-phase concentration", unit = u"mol/m^3"]
        c_s = 0.0, [description = "Surface concentration", unit = u"mol/m^3"]
    end

    @variables begin
        R_p(t), [description = "Particle radius", unit = u"m"]
    end

    eqs = [
        # Eq. 12.14: dR_p/dt = (D_g M_A)/(ρ_p R_p)(c_∞ - c_s)
        D(R_p) ~ (D_g * M_A / (ρ_p * R_p)) * (c_inf - c_s),
    ]

    return System(eqs, t; name)
end

"""
    MassTransferCoefficient(; name=:MassTransferCoefficient)

Calculate the combined gas-phase and interfacial mass transfer coefficient.

This coefficient accounts for both gas-phase diffusion and interfacial
resistance to mass transfer.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.116
"""
@component function MassTransferCoefficient(; name=:MassTransferCoefficient)
    @constants begin
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
    end

    @parameters begin
        R_p = 1.0e-5, [description = "Particle/droplet radius", unit = u"m"]
        D_g = 2.0e-5, [description = "Gas-phase diffusivity", unit = u"m^2/s"]
        α = 1.0, [description = "Accommodation coefficient (dimensionless)", unit = u"1"]
        M_A = 0.029, [description = "Molecular weight of species A", unit = u"kg/mol"]
        T = 298.15, [description = "Temperature", unit = u"K"]
    end

    @variables begin
        k_mt(t), [description = "Mass transfer coefficient", unit = u"s^-1"]
        τ_g(t), [description = "Gas-phase diffusion resistance", unit = u"s"]
        τ_i(t), [description = "Interfacial resistance", unit = u"s"]
    end

    eqs = [
        # Gas-phase diffusion resistance (first term in Eq. 12.116)
        τ_g ~ R_p^2 / (3 * D_g),
        # Interfacial resistance (second term in Eq. 12.116)
        τ_i ~ (R_p / (3 * α)) * sqrt(2 * π * M_A / (R * T)),
        # Eq. 12.116: k_mt = [R_p²/(3D_g) + (R_p/3α)√(2πM_A/(RT))]⁻¹
        k_mt ~ 1 / (τ_g + τ_i),
    ]

    return System(eqs, t; name)
end

"""
    UptakeCoefficient(; name=:UptakeCoefficient)

Calculate the uptake coefficient γ for gas-particle mass transfer.

The uptake coefficient is the probability that a molecule striking the
particle surface is taken up by the particle.

Reference: Seinfeld & Pandis (2006) Chapter 12, Eq. 12.121
"""
@component function UptakeCoefficient(; name=:UptakeCoefficient)
    @constants begin
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
    end

    @parameters begin
        R_p = 1.0e-5, [description = "Particle/droplet radius", unit = u"m"]
        D_g = 2.0e-5, [description = "Gas-phase diffusivity", unit = u"m^2/s"]
        α = 1.0, [description = "Accommodation coefficient (dimensionless)", unit = u"1"]
        c_bar = 300.0, [description = "Mean molecular speed", unit = u"m/s"]
        T = 298.15, [description = "Temperature", unit = u"K"]
        H_star = 1.0e5, [description = "Effective Henry's law coefficient", unit = u"mol/(m^3*Pa)"]
        k_rxn = 1.0, [description = "First-order reaction rate constant", unit = u"s^-1"]
        D_aq = 1.0e-9, [description = "Aqueous-phase diffusivity", unit = u"m^2/s"]
    end

    @variables begin
        γ(t), [description = "Uptake coefficient (dimensionless)", unit = u"1"]
        term_diff(t), [description = "Gas-phase diffusion term (dimensionless)", unit = u"1"]
        term_acc(t), [description = "Accommodation term (dimensionless)", unit = u"1"]
        term_rxn(t), [description = "Reaction term (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 12.121: 1/γ = (R_p c̄_A)/(4D_g) + 1/α + c̄_A/(4RT H*_A √(kD_{aq}))
        term_diff ~ (R_p * c_bar) / (4 * D_g),
        term_acc ~ 1 / α,
        term_rxn ~ c_bar / (4 * R * T * H_star * sqrt(k_rxn * D_aq)),
        γ ~ 1 / (term_diff + term_acc + term_rxn),
    ]

    return System(eqs, t; name)
end

"""
    MassTransfer(; name=:MassTransfer)

Complete mass transfer model combining transition regime correction,
gas-phase diffusion, and interfacial transport.

This component provides the transition-regime corrected mass flux
for particle growth/evaporation, applicable across all Knudsen numbers.

Reference: Seinfeld & Pandis (2006) Chapter 12, Sections 12.1.1-12.1.4
"""
@component function MassTransfer(; name=:MassTransfer)
    @constants begin
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
        k_B = 1.380649e-23, [description = "Boltzmann constant", unit = u"J/K"]
        N_Av = 6.02214076e23, [description = "Avogadro's number", unit = u"mol^-1"]
    end

    @parameters begin
        T = 298.15, [description = "Temperature", unit = u"K"]
        P = 101325.0, [description = "Pressure", unit = u"Pa"]
        M_A = 0.029, [description = "Molecular weight of species A", unit = u"kg/mol"]
        M_air = 0.029, [description = "Molecular weight of air", unit = u"kg/mol"]
        μ = 1.81e-5, [description = "Dynamic viscosity of air", unit = u"Pa*s"]
        R_p = 1.0e-7, [description = "Particle radius", unit = u"m"]
        D_g = 2.0e-5, [description = "Gas-phase diffusivity", unit = u"m^2/s"]
        α = 1.0, [description = "Accommodation coefficient (dimensionless)", unit = u"1"]
        c_inf = 1.0e-6, [description = "Bulk gas-phase concentration", unit = u"mol/m^3"]
        c_s = 0.0, [description = "Surface concentration", unit = u"mol/m^3"]
    end

    @variables begin
        c_bar(t), [description = "Mean molecular speed", unit = u"m/s"]
        λ(t), [description = "Mean free path", unit = u"m"]
        ρ_air(t), [description = "Air density", unit = u"kg/m^3"]
        Kn(t), [description = "Knudsen number (dimensionless)", unit = u"1"]
        f_FS(t), [description = "Fuchs-Sutugin correction factor (dimensionless)", unit = u"1"]
        J_c(t), [description = "Continuum regime molar flow rate", unit = u"mol/s"]
        J(t), [description = "Transition regime corrected molar flow rate", unit = u"mol/s"]
    end

    eqs = [
        # Mean molecular speed (Eq. 12.24)
        c_bar ~ sqrt(8 * k_B * T * N_Av / (π * M_A)),

        # Air density from ideal gas law
        ρ_air ~ P * M_air / (R * T),

        # Mean free path
        λ ~ 2 * μ / (ρ_air * sqrt(8 * R * T / (π * M_air))),

        # Knudsen number (Eq. 12.31)
        Kn ~ λ / R_p,

        # Fuchs-Sutugin correction (Eq. 12.43)
        f_FS ~ 0.75 * α * (1 + Kn) / (Kn^2 + Kn + 0.283 * Kn * α + 0.75 * α),

        # Continuum regime flux (Eq. 12.12)
        J_c ~ 4 * π * R_p * D_g * (c_inf - c_s),

        # Transition regime corrected flux
        J ~ f_FS * J_c,
    ]

    return System(eqs, t; name)
end
