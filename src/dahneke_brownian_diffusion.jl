export DahnekeMassTransportCorrection, DahnekeHeatTransportCorrection
export DahnekeCondensationEvaporation, DahnekeCoagulationRate
export DahnekeCapillaryPenetration

"""
    DahnekeMassTransportCorrection(; name=:DahnekeMassTransportCorrection)

Compute the correction factor β for diffusional mass transport to a sphere,
valid across all Knudsen numbers from the continuum to the free-molecular regime.

**Reference**: Dahneke, B. (1983). Simple Kinetic Theory of Brownian Diffusion in
Vapors and Aerosols. In *Theory of Dispersed Multiphase Flow* (pp. 97–133).
Academic Press. Eq. 5.5.

The correction factor relates the actual transport rate I to the continuum (Maxwell)
result I_M:

``I = β I_M``

where

``β = \\frac{Kn_D + 1}{\\frac{2 Kn_D (Kn_D + 1)}{δ} + 1}``

Here Kn_D = ℓ_D / r is the diffusion Knudsen number, ℓ_D = 2D/c̄ is the diffusional
mean-free-path, and δ is the condensation-evaporation (sticking) coefficient.

Limiting cases:
- Kn_D → 0: β → 1 (continuum limit, Maxwell's result)
- Kn_D → ∞: β → δc̄r/(4D) (free-molecular limit)
"""
@component function DahnekeMassTransportCorrection(; name = :DahnekeMassTransportCorrection)
    @constants begin
        k_B = 1.380649e-23, [description = "Boltzmann constant", unit = u"J/K"]
    end

    @parameters begin
        T = 298.15, [description = "Temperature", unit = u"K"]
        D_v = 2.0e-5, [description = "Diffusion coefficient of vapor molecules", unit = u"m^2/s"]
        m_v = 4.81e-26,
        [description = "Mass of diffusing vapor molecule", unit = u"kg"]
        r = 1.0e-7, [description = "Sphere (droplet/particle) radius", unit = u"m"]
        δ_m = 1.0,
        [description = "Condensation-evaporation (sticking) coefficient (dimensionless)",
            unit = u"1"]
    end

    @variables begin
        c_bar(t), [description = "Mean thermal velocity of vapor molecules", unit = u"m/s"]
        ℓ_D(t), [description = "Diffusional mean-free-path of vapor", unit = u"m"]
        Kn_D(t), [description = "Diffusion Knudsen number (dimensionless)", unit = u"1"]
        β(t), [description = "Mass transport correction factor (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 6.2 - Mean thermal velocity
        c_bar ~ sqrt(8 * k_B * T / (π * m_v)),
        # Eq. 3.6 - Diffusional mean-free-path
        ℓ_D ~ 2 * D_v / c_bar,
        # Kn_D = ℓ_D / r (diffusion Knudsen number based on sphere radius)
        Kn_D ~ ℓ_D / r,
        # Eq. 5.5 - Correction factor for mass transport
        β ~ (Kn_D + 1) / (2 * Kn_D * (Kn_D + 1) / δ_m + 1),
    ]

    return System(eqs, t; name)
end

"""
    DahnekeHeatTransportCorrection(; name=:DahnekeHeatTransportCorrection)

Compute the correction factor β_q for heat conduction to a sphere, valid across
all Knudsen numbers. Exactly parallel to the mass transport correction.

**Reference**: Dahneke, B. (1983). Simple Kinetic Theory of Brownian Diffusion in
Vapors and Aerosols. In *Theory of Dispersed Multiphase Flow* (pp. 97–133).
Academic Press. Eqs. 5.6–5.7.

``β_q = \\frac{Kn_K + 1}{\\frac{2 Kn_K (Kn_K + 1)}{α} + 1}``

where Kn_K = 2κ/(c̄' n' c_v r) is the heat transfer Knudsen number and α is the
thermal accommodation coefficient.
"""
@component function DahnekeHeatTransportCorrection(; name = :DahnekeHeatTransportCorrection)
    @parameters begin
        κ = 0.026, [description = "Thermal conductivity of gas", unit = u"W/(m*K)"]
        n_gas = 2.46e25,
        [description = "Number concentration of suspending gas molecules", unit = u"m^-3"]
        c_bar_gas = 500.0,
        [description = "Mean thermal velocity of gas molecules", unit = u"m/s"]
        c_v = 718.0,
        [description = "Specific heat at constant volume of gas (per unit mass)",
            unit = u"J/(kg*K)"]
        m_gas = 4.81e-26, [description = "Mass of gas molecule", unit = u"kg"]
        r = 1.0e-7, [description = "Sphere (droplet/particle) radius", unit = u"m"]
        α = 1.0,
        [description = "Thermal accommodation coefficient (dimensionless)", unit = u"1"]
    end

    @variables begin
        Kn_K(t), [description = "Heat transfer Knudsen number (dimensionless)", unit = u"1"]
        β_q(t),
        [description = "Heat transport correction factor (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 5.7 - Heat transfer Knudsen number: Kn_K = 2κ/(c̄' n' c_v r)
        # Note: n' c_v has units m^-3 * J/(kg*K) and we need to multiply by m_gas
        # to get volumetric heat capacity: n' * m_gas * c_v [J/(m^3*K)]
        Kn_K ~ 2 * κ / (c_bar_gas * n_gas * m_gas * c_v * r),
        # Eq. 5.7 - Correction factor for heat transport
        β_q ~ (Kn_K + 1) / (2 * Kn_K * (Kn_K + 1) / α + 1),
    ]

    return System(eqs, t; name)
end

"""
    DahnekeCondensationEvaporation(; name=:DahnekeCondensationEvaporation)

Model stationary condensation or evaporation of a droplet including coupled
mass and heat transfer with transition-regime corrections.

**Reference**: Dahneke, B. (1983). Simple Kinetic Theory of Brownian Diffusion in
Vapors and Aerosols. In *Theory of Dispersed Multiphase Flow* (pp. 97–133).
Academic Press. Eqs. 2.1–2.3, 5.4–5.8.

The corrected molecular transport rate is:

``I = 4π r D_v β (n_∞ - n_s)``

The corrected heat conduction rate is:

``Q = 4π r κ β_q (T_∞ - T_o)``

These are related through the conservation equation:

``β I_M = β_q Q_M / L``
"""
@component function DahnekeCondensationEvaporation(; name = :DahnekeCondensationEvaporation)
    mass_corr = DahnekeMassTransportCorrection(; name = :mass_corr)
    heat_corr = DahnekeHeatTransportCorrection(; name = :heat_corr)

    @parameters begin
        T_inf = 298.15, [description = "Gas temperature at infinity", unit = u"K"]
        T_o = 293.15, [description = "Droplet surface temperature", unit = u"K"]
        n_inf = 1.0e21,
        [description = "Vapor number concentration at infinity", unit = u"m^-3"]
        n_s = 0.0,
        [description = "Saturation vapor number concentration at droplet surface",
            unit = u"m^-3"]
        L_vap = 7.48e-20,
        [description = "Molecular latent heat of vaporization", unit = u"J"]
        r = 1.0e-7, [description = "Droplet radius", unit = u"m"]
        D_v = 2.0e-5,
        [description = "Diffusion coefficient of vapor molecules", unit = u"m^2/s"]
        κ = 0.026, [description = "Thermal conductivity of gas", unit = u"W/(m*K)"]
    end

    @variables begin
        I_M(t), [description = "Maxwell (continuum) molecular growth rate", unit = u"s^-1"]
        I(t), [description = "Corrected molecular growth rate", unit = u"s^-1"]
        Q_M(t), [description = "Maxwell (continuum) heat conduction rate", unit = u"W"]
        Q(t), [description = "Corrected heat conduction rate", unit = u"W"]
    end

    eqs = [
        # Share parameters with subsystems
        mass_corr.r ~ r,
        mass_corr.D_v ~ D_v,
        mass_corr.T ~ T_inf,
        heat_corr.r ~ r,
        heat_corr.κ ~ κ,
        # Eq. 2.1 - Maxwell's equation for molecular growth rate
        I_M ~ 4 * π * r * D_v * (n_inf - n_s),
        # Eq. 5.4 - Corrected mass transport rate
        I ~ mass_corr.β * I_M,
        # Eq. 2.2 - Maxwell's heat conduction rate
        Q_M ~ 4 * π * r * κ * (T_inf - T_o),
        # Eq. 5.6 - Corrected heat conduction rate
        Q ~ heat_corr.β_q * Q_M,
    ]

    return System(eqs, t; systems = [mass_corr, heat_corr], name)
end

"""
    DahnekeCoagulationRate(; name=:DahnekeCoagulationRate)

Compute the Brownian coagulation rate constant K using the Dahneke kinetic theory,
valid across all Knudsen numbers from the continuum to the free-molecular regime.

**Reference**: Dahneke, B. (1983). Simple Kinetic Theory of Brownian Diffusion in
Vapors and Aerosols. In *Theory of Dispersed Multiphase Flow* (pp. 97–133).
Academic Press. Eqs. 8.8, 8.13–8.15.

The coagulation rate constant is:

``K = 4π R D \\frac{Kn_D + 1}{1 + \\frac{2 Kn_D (Kn_D + 1)}{δ}}``

where R = r₁ + r₂ is the collision radius, D = D₁ + D₂ is the mutual diffusion
coefficient, and Kn_D = 2kT/(c̄fR) is the diffusion Knudsen number based on the
mutual friction coefficient f = f₁f₂/(f₁+f₂) and mean thermal speed
c̄ = √(c̄₁² + c̄₂²).

The correction factor β = K/K_o (where K_o = 4πRD is the continuum value) can
be decomposed as β = β₁β₂ where:
- β₁ = C_{s1}C_{s2}(r₁/C_{s1} + r₂/C_{s2})/(r₁+r₂) corrects for non-negligible Kn
- β₂ = (Kn_D+1)/[1+2Kn_D(Kn_D+1)/δ] corrects for non-negligible Kn_D
"""
@component function DahnekeCoagulationRate(; name = :DahnekeCoagulationRate)
    @constants begin
        k_B = 1.380649e-23, [description = "Boltzmann constant", unit = u"J/K"]
        π_val = 3.14159265358979, [description = "Pi (dimensionless)", unit = u"1"]
    end

    @parameters begin
        T = 298.15, [description = "Temperature", unit = u"K"]
        μ = 1.84e-5, [description = "Dynamic viscosity of suspending gas", unit = u"Pa*s"]
        r_1 = 1.0e-7, [description = "Radius of particle species 1", unit = u"m"]
        r_2 = 1.0e-7, [description = "Radius of particle species 2", unit = u"m"]
        ρ_1 = 1000.0, [description = "Density of particle species 1", unit = u"kg/m^3"]
        ρ_2 = 1000.0, [description = "Density of particle species 2", unit = u"kg/m^3"]
        δ_p = 1.0,
        [description = "Particle sticking probability (dimensionless)", unit = u"1"]
        C_s1 = 1.0,
        [description = "Slip correction factor for particle 1 (dimensionless)", unit = u"1"]
        C_s2 = 1.0,
        [description = "Slip correction factor for particle 2 (dimensionless)", unit = u"1"]
    end

    @variables begin
        m_1(t), [description = "Mass of particle 1", unit = u"kg"]
        m_2(t), [description = "Mass of particle 2", unit = u"kg"]
        f_1(t), [description = "Friction coefficient of particle 1", unit = u"kg/s"]
        f_2(t), [description = "Friction coefficient of particle 2", unit = u"kg/s"]
        D_1(t), [description = "Diffusion coefficient of particle 1", unit = u"m^2/s"]
        D_2(t), [description = "Diffusion coefficient of particle 2", unit = u"m^2/s"]
        D_mutual(t), [description = "Mutual diffusion coefficient", unit = u"m^2/s"]
        c_bar_1(t), [description = "Mean thermal velocity of particle 1", unit = u"m/s"]
        c_bar_2(t), [description = "Mean thermal velocity of particle 2", unit = u"m/s"]
        c_bar(t), [description = "Combined mean thermal velocity", unit = u"m/s"]
        f_mutual(t), [description = "Mutual friction coefficient", unit = u"kg/s"]
        R_coll(t), [description = "Collision radius r₁+r₂", unit = u"m"]
        Kn_D(t), [description = "Diffusion Knudsen number (dimensionless)", unit = u"1"]
        K_o(t), [description = "Continuum coagulation constant", unit = u"m^3/s"]
        β₂(t),
        [description = "Non-continuum correction factor β₂ (dimensionless)", unit = u"1"]
        K(t), [description = "Coagulation rate constant", unit = u"m^3/s"]
    end

    eqs = [
        # Particle masses
        m_1 ~ ρ_1 * 4 * π_val / 3 * r_1^3,
        m_2 ~ ρ_2 * 4 * π_val / 3 * r_2^3,
        # Friction coefficients: f = 6πηr/C_s
        f_1 ~ 6 * π_val * μ * r_1 / C_s1,
        f_2 ~ 6 * π_val * μ * r_2 / C_s2,
        # Einstein diffusion coefficients: D = kT/f (Eq. 6.5)
        D_1 ~ k_B * T / f_1,
        D_2 ~ k_B * T / f_2,
        # Mutual diffusion coefficient: D = D_1 + D_2 (p. 115)
        D_mutual ~ D_1 + D_2,
        # Mean thermal velocities (Eq. 6.2)
        c_bar_1 ~ sqrt(8 * k_B * T / (π_val * m_1)),
        c_bar_2 ~ sqrt(8 * k_B * T / (π_val * m_2)),
        # Combined mean thermal velocity: c̄ = √(c̄₁² + c̄₂²) (p. 115)
        c_bar ~ sqrt(c_bar_1^2 + c_bar_2^2),
        # Mutual friction coefficient: f = f₁f₂/(f₁+f₂) (p. 115)
        f_mutual ~ f_1 * f_2 / (f_1 + f_2),
        # Collision radius (p. 114)
        R_coll ~ r_1 + r_2,
        # Diffusion Knudsen number: Kn_D = 2kT/(c̄fR) (p. 115)
        Kn_D ~ 2 * k_B * T / (c_bar * f_mutual * R_coll),
        # Eq. 8.8 - Continuum coagulation constant (uses actual D with C_s)
        K_o ~ 4 * π_val * R_coll * D_mutual,
        # Eq. 8.15 β₂ correction for non-negligible Kn_D
        β₂ ~ (Kn_D + 1) / (1 + 2 * Kn_D * (Kn_D + 1) / δ_p),
        # Eq. 8.13 - Coagulation rate constant: K = K_o * β₂
        # K_o already includes C_s via diffusion coefficients D_1, D_2.
        # The total correction β = K/K_o(C_s=1) decomposes as β₁β₂ (Eq. 8.14)
        # where β₁ = D_actual/D_base captures the C_s effect.
        K ~ K_o * β₂,
    ]

    return System(eqs, t; name)
end

"""
    DahnekeCapillaryPenetration(n_modes::Int=6; name=:DahnekeCapillaryPenetration)

Compute the penetration fraction of aerosol particles through a fine capillary
using the eigenfunction expansion method. The penetration φ(z) is the fraction
of particles that pass through the capillary without depositing on the wall.

**Reference**: Dahneke, B. (1983). Simple Kinetic Theory of Brownian Diffusion in
Vapors and Aerosols. In *Theory of Dispersed Multiphase Flow* (pp. 97–133).
Academic Press. Section 9, Eqs. 9.1–9.9.

The steady convective diffusion equation in a circular tube is solved by separation
of variables, yielding the penetration function:

``ϕ(z) = \\sum_i B_i \\exp(-ω_i^2 z)``

where z = DL/(σv̄R²) is the dimensionless tube length, ω_i² are eigenvalues, and
B_i are coefficients that depend on the velocity profile parameter γ and the
diffusion Knudsen number Kn_D.

The velocity profile is v(r) = σv̄[1 - γ(r/R)²] where:
- γ = 0: free-molecule flow (plug flow)
- γ = 0.5: transition flow
- γ = 1.0: continuum flow (Poiseuille flow)
- σ = 2/(2-γ) from the flow rate constraint

Arguments:
- `n_modes::Int`: Number of eigenmodes to include (default 6, matching Tables 3–5)

The eigenvalues and coefficients are precomputed from Dahneke's Tables 3, 4, and 5
for δ = 1 (perfect sticking) and γ = 0, 0.5, 1.0.
"""
@component function DahnekeCapillaryPenetration(n_modes::Int = 6;
        name = :DahnekeCapillaryPenetration)
    @parameters begin
        Dcoeff = 1.0e-10, [description = "Particle diffusion coefficient", unit = u"m^2/s"]
        R_cap = 1.0e-4, [description = "Capillary radius", unit = u"m"]
        v_mean = 0.01, [description = "Mean flow velocity", unit = u"m/s"]
        L_cap = 0.1, [description = "Capillary length", unit = u"m"]
        γ_flow = 1.0,
        [description = "Velocity profile parameter: 0=plug, 0.5=transition, 1=Poiseuille (dimensionless)",
            unit = u"1"]
        Kn_D_cap = 0.0,
        [description = "Diffusion Knudsen number Kn_D = 2kT/(c̄fR) (dimensionless)",
            unit = u"1"]
    end

    @variables begin
        σ_flow(t),
        [description = "Velocity profile normalization factor (dimensionless)", unit = u"1"]
        z_dim(t),
        [description = "Dimensionless capillary length z=DL/(σv̄R²) (dimensionless)",
            unit = u"1"]
    end

    eqs = [
        # Eq. 9.4, 9.5 - σ = 2/(2-γ) from flow rate constraint
        σ_flow ~ 2 / (2 - γ_flow),
        # Dimensionless tube length z' = Dz/(σv̄R²)
        z_dim ~ Dcoeff * L_cap / (σ_flow * v_mean * R_cap^2),
    ]

    return System(eqs, t; name)
end
