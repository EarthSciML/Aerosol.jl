export DahnekeMassTransportCorrection, DahnekeHeatTransportCorrection
export DahnekeCondensationEvaporation, DahnekeCoagulationRate
export DahnekeCapillaryPenetration, dahneke_capillary_penetration
export DAHNEKE_TABLE3, DAHNEKE_TABLE4, DAHNEKE_TABLE5

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
        D_v = 2.0e-5,
        [description = "Diffusion coefficient of vapor molecules", unit = u"m^2/s"]
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
        β(t),
        [description = "Mass transport correction factor (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 6.2 - Mean thermal velocity
        c_bar ~ sqrt(8 * k_B * T / (π * m_v)),
        # Eq. 3.6 - Diffusional mean-free-path
        ℓ_D ~ 2 * D_v / c_bar,
        # Kn_D = ℓ_D / r (diffusion Knudsen number based on sphere radius)
        Kn_D ~ ℓ_D / r,
        # Eq. 5.5 - Correction factor for mass transport
        β ~ (Kn_D + 1) / (2 * Kn_D * (Kn_D + 1) / δ_m + 1)
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
        β_q ~ (Kn_K + 1) / (2 * Kn_K * (Kn_K + 1) / α + 1)
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
        Q ~ heat_corr.β_q * Q_M
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
        K ~ K_o * β₂
    ]

    return System(eqs, t; name)
end

"""
    DahnekeCapillaryPenetration(; name=:DahnekeCapillaryPenetration)

Compute parameters for aerosol penetration through a fine capillary tube.

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
  - σ = 2/(2-γ) from the flow rate constraint (Eqs. 9.4–9.5)

The eigenvalues ω_i² and coefficients B_i are tabulated in Tables 3 (γ=0),
4 (γ=0.5), and 5 (γ=1.0) of the paper for δ=1 and various Kn_D values.
Use [`dahneke_capillary_penetration`](@ref) to compute the penetration fraction
from these tables.
"""
@component function DahnekeCapillaryPenetration(; name = :DahnekeCapillaryPenetration)
    @parameters begin
        Dcoeff = 1.0e-10, [description = "Particle diffusion coefficient", unit = u"m^2/s"]
        R_cap = 1.0e-4, [description = "Capillary radius", unit = u"m"]
        v_mean = 0.01, [description = "Mean flow velocity", unit = u"m/s"]
        L_cap = 0.1, [description = "Capillary length", unit = u"m"]
        γ_flow = 1.0,
        [
            description = "Velocity profile parameter: 0=plug, 0.5=transition, 1=Poiseuille (dimensionless)",
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
        z_dim ~ Dcoeff * L_cap / (σ_flow * v_mean * R_cap^2)
    ]

    return System(eqs, t; name)
end

# Eigenvalue tables from Dahneke (1983) Tables 3, 4, and 5.
# Each table maps Kn_D => (ω_i², B_i) for i=1..6 eigenmodes with δ=1.

"""
Eigenvalues ω_i² and coefficients B_i from Table 3 of Dahneke (1983)
for γ = 0 (free-molecule / plug flow) with δ = 1.
Keys are Kn_D values; values are vectors of (ω_i², B_i) tuples for i = 1..6.
"""
const DAHNEKE_TABLE3 = Dict(
    0.0 => [
        (5.78319, 0.69166), (30.4713, 0.13127), (74.8870, 0.05341),
        (139.040, 0.02877), (222.932, 0.01794), (326.563, 0.01125)],
    0.1 => [
        (4.75021, 0.80388), (25.3332, 0.12598), (63.3120, 0.03869),
        (119.603, 0.01523), (194.827, 0.00696), (289.336, 0.00353)],
    0.2 => [
        (3.95936, 0.87214), (22.2137, 0.09535), (58.0295, 0.02075),
        (112.833, 0.00642), (187.103, 0.00252), (280.998, 0.00116)],
    0.3 => [
        (3.36390, 0.91276), (20.3669, 0.06932), (55.4762, 0.01203),
        (109.951, 0.00334), (184.057, 0.00124), (277.861, 0.00055)],
    0.4 => [
        (2.91051, 0.93767), (19.2003, 0.05116), (54.0341, 0.00767),
        (108.406, 0.00201), (182.464, 0.00073), (276.244, 0.00032)],
    0.5 => [
        (2.55824, 0.95366), (18.4123, 0.03877), (53.1206, 0.00527),
        (107.450, 0.00134), (181.492, 0.00048), (275.262, 0.00021)]
)

"""
Eigenvalues ω_i² and coefficients B_i from Table 4 of Dahneke (1983)
for γ = 0.5 (transition / slip flow) with δ = 1.
Keys are Kn_D values; values are vectors of (ω_i², B_i) tuples for i = 1..6.
"""
const DAHNEKE_TABLE4 = Dict(
    0.0 => [
        (6.47641, 0.72680), (36.1924, 0.11940), (89.8982, 0.04708),
        (167.527, 0.02505), (269.062, 0.01536), (394.496, 0.00020)],
    0.1 => [
        (5.45150, 0.81527), (31.1557, 0.11352), (78.3677, 0.03671),
        (147.702, 0.01558), (239.731, 0.00758), (354.873, 0.00098)],
    0.2 => [
        (4.64778, 0.87221), (27.8208, 0.09105), (72.2133, 0.02230),
        (139.209, 0.00757), (229.500, 0.00314), (343.394, 0.00039)],
    0.3 => [
        (4.02261, 0.90850), (25.6742, 0.07000), (68.9119, 0.01392),
        (135.226, 0.00418), (225.124, 0.00161), (338.782, 0.00020)],
    0.4 => [
        (3.53186, 0.93215), (24.2373, 0.05397), (66.9437, 0.00927),
        (133.002, 0.00260), (222.769, 0.00097), (336.354, 0.00013)],
    0.5 => [
        (3.14060, 0.94809), (23.2277, 0.04230), (65.6588, 0.00655),
        (131.600, 0.00176), (221.310, 0.00064), (334.864, 0.00008)]
)

"""
Eigenvalues ω_i² and coefficients B_i from Table 5 of Dahneke (1983)
for γ = 1.0 (continuum / Poiseuille flow) with δ = 1.
Keys are Kn_D values; values are vectors of (ω_i², B_i) tuples for i = 1..6.
"""
const DAHNEKE_TABLE5 = Dict(
    0.0 => [
        (7.31359, 0.81905), (44.6095, 0.09753), (113.921, 0.03250),
        (215.241, 0.01543), (348.564, 0.00168), (513.890, 0.00000)],
    0.1 => [
        (6.33404, 0.86923), (40.5081, 0.08320), (105.487, 0.02358),
        (201.646, 0.00977), (329.192, 0.00379), (488.261, 0.00000)],
    0.2 => [
        (5.55381, 0.90293), (37.6387, 0.06760), (100.270, 0.01633),
        (194.094, 0.00602), (319.401, 0.00217), (476.348, 0.00000)],
    0.3 => [
        (4.92768, 0.92594), (35.6031, 0.05441), (96.9237, 0.01156),
        (189.606, 0.00393), (313.917, 0.00135), (469.977, 0.00000)],
    0.4 => [
        (4.41898, 0.94205), (34.1151, 0.04405), (94.6517, 0.00847),
        (186.705, 0.00273), (310.492, 0.00091), (466.107, 0.00000)],
    0.5 => [
        (4.00000, 0.95362), (32.9926, 0.03608), (93.0272, 0.00642),
        (184.697, 0.00199), (308.171, 0.00066), (463.521, 0.00000)]
)

"""
    dahneke_capillary_penetration(z, γ, Kn_D)

Compute the aerosol penetration fraction ϕ through a capillary tube at
dimensionless length z, for velocity profile parameter γ and diffusion
Knudsen number Kn_D, using the eigenfunction expansion from Dahneke (1983)
Section 9.

The penetration is computed as:

``ϕ(z) = \\sum_{i=1}^{6} B_i \\exp(-ω_i^2 z)``

where the eigenvalues ω_i² and coefficients B_i are taken from Tables 3–5
of the paper for δ = 1 (perfect sticking).

# Arguments

  - `z`: Dimensionless tube length z = DL/(σv̄R²)
  - `γ`: Velocity profile parameter (0=plug, 0.5=transition, 1.0=Poiseuille)
  - `Kn_D`: Diffusion Knudsen number (must be one of 0.0, 0.1, 0.2, 0.3, 0.4, 0.5)

# Returns

  - Penetration fraction ϕ (between 0 and 1)
"""
function dahneke_capillary_penetration(z, γ, Kn_D)
    if γ ≈ 0.0
        table = DAHNEKE_TABLE3
    elseif γ ≈ 0.5
        table = DAHNEKE_TABLE4
    elseif γ ≈ 1.0
        table = DAHNEKE_TABLE5
    else
        error("γ must be 0.0, 0.5, or 1.0; got γ = $γ")
    end

    if !haskey(table, Kn_D)
        error("Kn_D must be one of $(sort(collect(keys(table)))); got Kn_D = $Kn_D")
    end

    modes = table[Kn_D]
    ϕ = sum(B * exp(-ω² * z) for (ω², B) in modes)
    return ϕ
end
