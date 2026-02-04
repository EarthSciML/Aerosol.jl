export AerosolLayerRadiativeForcing, CriticalSingleScatteringAlbedo,
       CloudOpticalDepth, CloudAlbedo, CloudAlbedoSensitivity, IndirectAerosolForcing

"""
    AerosolLayerRadiativeForcing(; name=:AerosolLayerRadiativeForcing)

Aerosol layer radiative forcing model from Seinfeld & Pandis (2006) Chapter 24.

This component calculates the change in outgoing radiative flux due to an aerosol
layer using the scattering-absorbing model from Section 24.1. It computes the
reflection and transmission coefficients, the total reflectance of the aerosol-surface
system, and the change in planetary albedo.

**Reference**: Seinfeld, J. H., and Pandis, S. N. (2006). *Atmospheric Chemistry and
Physics: From Air Pollution to Climate Change*, 2nd Edition, Chapter 24, Equations
24.1-24.10.

# Parameters

  - `F_0`: Incident solar flux (W m⁻²)
  - `τ`: Aerosol optical depth (dimensionless)
  - `ω`: Single-scattering albedo (dimensionless)
  - `β`: Upscatter fraction (dimensionless)
  - `R_s`: Surface albedo (dimensionless)
  - `A_c`: Cloud fraction (dimensionless)
  - `T_a`: Atmospheric transmittance (dimensionless)

# Variables (computed)

  - `r_aer`: Fraction of light reflected upward by aerosol (dimensionless)
  - `t_aer`: Fraction of light transmitted through aerosol (dimensionless)
  - `R_as`: Total reflectance of aerosol-surface system (dimensionless)
  - `ΔR_p`: Change in planetary albedo (dimensionless)
  - `ΔF`: Change in outgoing radiative flux (W m⁻²)
"""
@component function AerosolLayerRadiativeForcing(; name = :AerosolLayerRadiativeForcing)
    @parameters begin
        F_0 = 1370.0, [unit = u"W/m^2", description = "Incident solar flux"]
        τ = 0.1, [unit = u"1", description = "Aerosol optical depth (dimensionless)"]
        ω = 0.95, [unit = u"1", description = "Single-scattering albedo (dimensionless)"]
        β = 0.29, [unit = u"1", description = "Upscatter fraction (dimensionless)"]
        R_s = 0.15, [unit = u"1", description = "Surface albedo (dimensionless)"]
        A_c = 0.6, [unit = u"1", description = "Cloud fraction (dimensionless)"]
        T_a = 0.76, [unit = u"1", description = "Atmospheric transmittance (dimensionless)"]
    end

    @variables begin
        r_aer(t),
        [unit = u"1", description = "Fraction of light reflected upward by aerosol"]
        t_aer(t),
        [unit = u"1", description = "Fraction of light transmitted through aerosol"]
        R_as(t), [unit = u"1", description = "Total reflectance of aerosol-surface system"]
        ΔR_p(t), [unit = u"1", description = "Change in planetary albedo"]
        ΔF(t), [unit = u"W/m^2", description = "Change in outgoing radiative flux"]
    end

    eqs = [
        # Eq. 24.1 auxiliary: Fraction reflected upward by aerosol layer
        # r = (1 - exp(-τ)) * ω * β
        r_aer ~ (1 - exp(-τ)) * ω * β,

        # Eq. 24.1: Total fraction transmitted downward
        # t = exp(-τ) + ω*(1-β)*(1-exp(-τ))
        t_aer ~ exp(-τ) + ω * (1 - β) * (1 - exp(-τ)),

        # Eq. 24.5: Total reflectance of aerosol-surface system
        # R_as = r + t²*R_s / (1 - R_s*r)
        R_as ~ r_aer + t_aer^2 * R_s / (1 - R_s * r_aer),

        # Eq. 24.9: Change in planetary albedo (accounting for clouds and atmosphere)
        # ΔR_p = (1 - A_c) * T_a² * (R_as - R_s)
        ΔR_p ~ (1 - A_c) * T_a^2 * (R_as - R_s),

        # Eq. 24.10: Change in outgoing radiative flux
        # ΔF = F_0 * ΔR_p
        ΔF ~ F_0 * ΔR_p
    ]

    return System(eqs, t; name)
end

"""
    CriticalSingleScatteringAlbedo(; name=:CriticalSingleScatteringAlbedo)

Component to calculate the critical single-scattering albedo (ω_crit) at which an
aerosol layer transitions from cooling to heating the climate system.

From Seinfeld & Pandis (2006) Chapter 24, Equation 24.15:
ω_crit = 2*R_s / (2*R_s + β*(1-R_s)²)

When ω > ω_crit, the aerosol causes net cooling (increased reflection).
When ω < ω_crit, the aerosol causes net heating (absorption dominates).

**Reference**: Seinfeld, J. H., and Pandis, S. N. (2006). *Atmospheric Chemistry and
Physics*, Chapter 24, Section 24.2.

# Parameters

  - `R_s`: Surface albedo (dimensionless)
  - `β`: Upscatter fraction (dimensionless)

# Variables (computed)

  - `ω_crit`: Critical single-scattering albedo (dimensionless)
"""
@component function CriticalSingleScatteringAlbedo(; name = :CriticalSingleScatteringAlbedo)
    @parameters begin
        R_s = 0.15, [unit = u"1", description = "Surface albedo (dimensionless)"]
        β = 0.29, [unit = u"1", description = "Upscatter fraction (dimensionless)"]
    end

    @variables begin
        ω_crit(t),
        [unit = u"1", description = "Critical single-scattering albedo (dimensionless)"]
    end

    eqs = [
    # Eq. 24.15: Critical single-scattering albedo
        ω_crit ~ 2 * R_s / (2 * R_s + β * (1 - R_s)^2)
    ]

    return System(eqs, t; name)
end

"""
    CloudOpticalDepth(; name=:CloudOpticalDepth)

Component to calculate cloud optical depth from liquid water content and droplet
number concentration.

From Seinfeld & Pandis (2006) Chapter 24, Section 24.8.1:
τ_c = h * (9π*L²*N / 2ρ_w²)^(1/3)    (Eq. 24.36)

This relates cloud optical depth to cloud thickness (h), liquid water content (L),
and cloud droplet number concentration (N).

**Reference**: Seinfeld, J. H., and Pandis, S. N. (2006). *Atmospheric Chemistry and
Physics*, Chapter 24, Section 24.8.1, Equations 24.34-24.36.

# Parameters

  - `h`: Cloud thickness (m)
  - `L`: Liquid water content (kg m⁻³)
  - `N`: Cloud droplet number concentration (m⁻³)

# Constants

  - `ρ_w`: Density of water (kg m⁻³)

# Variables (computed)

  - `τ_c`: Cloud optical depth (dimensionless)
"""
@component function CloudOpticalDepth(; name = :CloudOpticalDepth)
    @constants begin
        ρ_w = 1000.0, [unit = u"kg/m^3", description = "Density of water"]
    end

    @parameters begin
        h = 500.0, [unit = u"m", description = "Cloud thickness"]
        L = 0.3e-3, [unit = u"kg/m^3", description = "Liquid water content"]
        N = 100.0e6, [unit = u"m^-3", description = "Cloud droplet number concentration"]
    end

    @variables begin
        τ_c(t), [unit = u"1", description = "Cloud optical depth (dimensionless)"]
    end

    eqs = [
    # Eq. 24.36: Cloud optical depth in terms of L, h, N
    # τ_c = h * (9π*L²*N / 2ρ_w²)^(1/3)
        τ_c ~ h * (9 * π * L^2 * N / (2 * ρ_w^2))^(1/3)
    ]

    return System(eqs, t; name)
end

"""
    CloudAlbedo(; name=:CloudAlbedo)

Component to calculate cloud albedo from optical depth using the two-stream
approximation.

From Seinfeld & Pandis (2006) Chapter 24, Section 24.8.1:
R_c = τ_c / (τ_c + γ)    (Eq. 24.38, simplified form for g = 0.85)

where γ = 7.7 for an asymmetry factor g = 0.85 (typical for cloud droplets).

The general two-stream approximation (Eq. 24.37) is:
R_c = √3*(1-g)*τ_c / (2 + √3*(1-g)*τ_c)

**Reference**: Seinfeld, J. H., and Pandis, S. N. (2006). *Atmospheric Chemistry and
Physics*, Chapter 24, Section 24.8.1, Equations 24.37-24.38.

# Parameters

  - `τ_c`: Cloud optical depth (dimensionless)
  - `g`: Asymmetry factor (dimensionless, default 0.85)

# Variables (computed)

  - `R_c`: Cloud albedo (dimensionless)
  - `γ`: Denominator coefficient for two-stream approximation (dimensionless)
"""
@component function CloudAlbedo(; name = :CloudAlbedo)
    @parameters begin
        τ_c = 10.0, [unit = u"1", description = "Cloud optical depth (dimensionless)"]
        g = 0.85, [unit = u"1", description = "Asymmetry factor (dimensionless)"]
    end

    @variables begin
        R_c(t), [unit = u"1", description = "Cloud albedo (dimensionless)"]
        γ(t),
        [unit = u"1", description = "Two-stream approximation coefficient (dimensionless)"]
    end

    eqs = [
        # Coefficient in two-stream approximation: γ = 2/(√3*(1-g))
        # For g = 0.85, γ ≈ 7.7
        γ ~ 2 / (sqrt(3) * (1 - g)),

        # Eq. 24.37/24.38: Cloud albedo (two-stream approximation)
        R_c ~ τ_c / (τ_c + γ)
    ]

    return System(eqs, t; name)
end

"""
    CloudAlbedoSensitivity(; name=:CloudAlbedoSensitivity)

Component to calculate the sensitivity of cloud albedo to cloud droplet number
concentration (Twomey susceptibility).

From Seinfeld & Pandis (2006) Chapter 24, Section 24.8.2:
dR_c/dN = R_c*(1-R_c)/(3N)    (Eq. 24.40)
dR_c/d(ln N) = R_c*(1-R_c)/3  (Eq. 24.41, Twomey susceptibility)

The maximum susceptibility occurs at R_c = 0.5.

**Reference**: Seinfeld, J. H., and Pandis, S. N. (2006). *Atmospheric Chemistry and
Physics*, Chapter 24, Section 24.8.2, Equations 24.39-24.42.

# Parameters

  - `R_c`: Cloud albedo (dimensionless)
  - `N`: Cloud droplet number concentration (m⁻³)

# Variables (computed)

  - `dRc_dN`: Sensitivity of cloud albedo to CDNC (m³)
  - `S`: Twomey susceptibility dR_c/d(ln N) (dimensionless)
"""
@component function CloudAlbedoSensitivity(; name = :CloudAlbedoSensitivity)
    @parameters begin
        R_c = 0.5, [unit = u"1", description = "Cloud albedo (dimensionless)"]
        N = 100.0e6, [unit = u"m^-3", description = "Cloud droplet number concentration"]
    end

    @variables begin
        dRc_dN(t), [unit = u"m^3", description = "Sensitivity of cloud albedo to CDNC"]
        S(t),
        [unit = u"1", description = "Twomey susceptibility dR_c/d(ln N) (dimensionless)"]
    end

    eqs = [
        # Eq. 24.40: Cloud albedo sensitivity to N
        dRc_dN ~ R_c * (1 - R_c) / (3 * N),

        # Eq. 24.41: Twomey susceptibility
        S ~ R_c * (1 - R_c) / 3
    ]

    return System(eqs, t; name)
end

"""
    IndirectAerosolForcing(; name=:IndirectAerosolForcing)

Component to calculate the indirect radiative forcing from aerosol effects on cloud
albedo (first indirect effect, or Twomey effect).

From Seinfeld & Pandis (2006) Chapter 24, Section 24.8.2:
ΔF_c = -F_0 * A_c * T_a² * ΔR_c    (Eq. 24.43)
ΔR_c = R_c*(1-R_c)/3 * Δln(N)      (Eq. 24.42)

This gives the change in shortwave forcing due to changes in cloud albedo caused
by changes in cloud droplet number concentration.

**Reference**: Seinfeld, J. H., and Pandis, S. N. (2006). *Atmospheric Chemistry and
Physics*, Chapter 24, Section 24.8.2, Equations 24.42-24.43.

# Parameters

  - `F_0`: Incident solar flux (W m⁻²)
  - `A_c`: Cloud fraction (dimensionless)
  - `T_a`: Atmospheric transmittance (dimensionless)
  - `R_c`: Cloud albedo (dimensionless)
  - `Δln_N`: Fractional change in CDNC, Δln(N) = ΔN/N (dimensionless)

# Variables (computed)

  - `ΔR_c`: Change in cloud albedo (dimensionless)
  - `ΔF_c`: Change in shortwave forcing (W m⁻²)
"""
@component function IndirectAerosolForcing(; name = :IndirectAerosolForcing)
    @parameters begin
        F_0 = 1370.0, [unit = u"W/m^2", description = "Incident solar flux"]
        A_c = 0.6, [unit = u"1", description = "Cloud fraction (dimensionless)"]
        T_a = 0.76, [unit = u"1", description = "Atmospheric transmittance (dimensionless)"]
        R_c = 0.5, [unit = u"1", description = "Cloud albedo (dimensionless)"]
        Δln_N = 0.262,
        [unit = u"1", description = "Fractional change in CDNC ln(N₁/N₀) (dimensionless)"]
    end

    @variables begin
        ΔR_c(t), [unit = u"1", description = "Change in cloud albedo (dimensionless)"]
        ΔF_c(t), [unit = u"W/m^2", description = "Change in shortwave forcing"]
    end

    eqs = [
        # Eq. 24.42: Change in cloud albedo
        ΔR_c ~ R_c * (1 - R_c) / 3 * Δln_N,

        # Eq. 24.43: Change in shortwave forcing (negative = cooling)
        ΔF_c ~ -F_0 * A_c * T_a^2 * ΔR_c
    ]

    return System(eqs, t; name)
end
