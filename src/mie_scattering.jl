export MieScattering, RayleighScattering, AerosolExtinction, Visibility

"""
    MieScattering(; name=:MieScattering)

Mie theory for scattering and absorption of light by spherical particles.

This implements the equations from Chapter 15 of Seinfeld & Pandis (2006)
"Atmospheric Chemistry and Physics: From Air Pollution to Climate Change".

The Mie theory calculates the scattering and extinction efficiencies (Q_scat, Q_ext)
as functions of the size parameter α = πDₚ/λ and the complex refractive index m = n + ik.

# Key Variables
- `Q_scat`: Scattering efficiency (dimensionless)
- `Q_abs`: Absorption efficiency (dimensionless)
- `Q_ext`: Extinction efficiency (dimensionless)
- `ω`: Single-scattering albedo (dimensionless)
- `g`: Asymmetry parameter (dimensionless)

# Key Parameters
- `D_p`: Particle diameter (m)
- `λ`: Wavelength of incident radiation (m)
- `n_refr`: Real part of refractive index (dimensionless)
- `k_refr`: Imaginary part of refractive index (dimensionless)
- `ρ_p`: Particle density (kg/m³)

# Reference
Seinfeld, J. H., & Pandis, S. N. (2006). Atmospheric Chemistry and Physics (2nd ed.),
Chapter 15, Equations 15.1-15.23, 15.37-15.43.
"""
@component function MieScattering(; name=:MieScattering)
    @constants begin
        π_const = Float64(π), [description = "Mathematical constant pi (dimensionless)"]
        length_scale = 1.0, [unit = u"m", description = "Length scale for dimensionless conversion"]
    end

    @parameters begin
        D_p = 0.5e-6, [unit = u"m", description = "Particle diameter"]
        λ = 550e-9, [unit = u"m", description = "Wavelength of incident radiation"]
        n_refr = 1.5, [description = "Real part of complex refractive index (dimensionless)"]
        k_refr = 0.0, [description = "Imaginary part of complex refractive index (dimensionless)"]
        ρ_p = 1500.0, [unit = u"kg/m^3", description = "Particle density"]
    end

    @variables begin
        α(t), [description = "Size parameter (dimensionless)"]
        Q_scat(t), [description = "Scattering efficiency (dimensionless)"]
        Q_abs(t), [description = "Absorption efficiency (dimensionless)"]
        Q_ext(t), [description = "Extinction efficiency (dimensionless)"]
        ω(t), [description = "Single-scattering albedo (dimensionless)"]
        C_scat(t), [unit = u"m^2", description = "Scattering cross section"]
        C_abs(t), [unit = u"m^2", description = "Absorption cross section"]
        C_ext(t), [unit = u"m^2", description = "Extinction cross section"]
        E_scat(t), [unit = u"m^2/kg", description = "Mass scattering efficiency"]
        E_abs(t), [unit = u"m^2/kg", description = "Mass absorption efficiency"]
        E_ext(t), [unit = u"m^2/kg", description = "Mass extinction efficiency"]
    end

    eqs = [
        # Eq. 15.6 - Size parameter (dimensionless ratio)
        α ~ π_const * (D_p / length_scale) / (λ / length_scale),

        # Mie efficiencies computed via registered functions (see below)
        Q_ext ~ mie_Q_ext(α, n_refr, k_refr),  # Eq. 15.14
        Q_scat ~ mie_Q_scat(α, n_refr, k_refr), # Eq. 15.13

        # Eq. 15.4 - Absorption efficiency from extinction and scattering
        Q_abs ~ Q_ext - Q_scat,

        # Eq. 15.5 - Single-scattering albedo
        ω ~ Q_scat / Q_ext,

        # Cross sections from efficiencies (C = Q × A, where A = πD_p²/4)
        C_scat ~ Q_scat * π_const * D_p^2 / 4,
        C_abs ~ Q_abs * π_const * D_p^2 / 4,
        C_ext ~ Q_ext * π_const * D_p^2 / 4,  # Eq. 15.3 implied

        # Eq. 15.41, 15.42, 15.43 - Mass efficiencies
        E_ext ~ 3 * Q_ext / (2 * ρ_p * D_p),
        E_scat ~ 3 * Q_scat / (2 * ρ_p * D_p),
        E_abs ~ 3 * Q_abs / (2 * ρ_p * D_p),
    ]

    return System(eqs, t; name, checks=false)
end

"""
    RayleighScattering(; name=:RayleighScattering)

Rayleigh scattering regime for particles much smaller than the wavelength (α << 1).

This provides closed-form solutions valid when D_p << λ (typically D_p ≤ 0.1 μm
for visible light). In this regime:
- Scattered intensity is proportional to 1/λ⁴
- Scattering pattern is symmetric in forward/backward directions
- Absorption dominates over scattering (Q_abs >> Q_scat)

# Reference
Seinfeld & Pandis (2006), Chapter 15, Equations 15.17-15.21.
"""
@component function RayleighScattering(; name=:RayleighScattering)
    @constants begin
        π_const = Float64(π), [description = "Mathematical constant pi (dimensionless)"]
        length_scale = 1.0, [unit = u"m", description = "Length scale for dimensionless conversion"]
    end

    @parameters begin
        D_p = 0.05e-6, [unit = u"m", description = "Particle diameter"]
        λ = 550e-9, [unit = u"m", description = "Wavelength of incident radiation"]
        n_refr = 1.5, [description = "Real part of complex refractive index (dimensionless)"]
        k_refr = 0.0, [description = "Imaginary part of complex refractive index (dimensionless)"]
        ρ_p = 1500.0, [unit = u"kg/m^3", description = "Particle density"]
    end

    @variables begin
        α(t), [description = "Size parameter (dimensionless)"]
        m2_factor_real(t), [description = "Real part of (m²-1)/(m²+2) factor (dimensionless)"]
        m2_factor_imag(t), [description = "Imaginary part of (m²-1)/(m²+2) factor (dimensionless)"]
        m2_factor_mag_sq(t), [description = "|(m²-1)/(m²+2)|² factor (dimensionless)"]
        Q_scat(t), [description = "Rayleigh scattering efficiency (dimensionless)"]
        Q_abs(t), [description = "Rayleigh absorption efficiency (dimensionless)"]
        Q_ext(t), [description = "Rayleigh extinction efficiency (dimensionless)"]
        ω(t), [description = "Single-scattering albedo (dimensionless)"]
        E_scat(t), [unit = u"m^2/kg", description = "Mass scattering efficiency"]
        E_abs(t), [unit = u"m^2/kg", description = "Mass absorption efficiency"]
        E_ext(t), [unit = u"m^2/kg", description = "Mass extinction efficiency"]
    end

    eqs = [
        # Size parameter (Eq. 15.6) - dimensionless ratio
        α ~ π_const * (D_p / length_scale) / (λ / length_scale),

        # Compute the (m²-1)/(m²+2) factor for complex m = n + ik
        # m² = n² - k² + 2nik
        # Let A = n² - k² - 1, B = 2nk, C = n² - k² + 2, D = 2nk
        # Then (m²-1)/(m²+2) = (A + iB)/(C + iD)
        # Real part = (AC + BD)/(C² + D²)
        # Imag part = (BC - AD)/(C² + D²)
        m2_factor_real ~ rayleigh_m2_factor_real(n_refr, k_refr),
        m2_factor_imag ~ rayleigh_m2_factor_imag(n_refr, k_refr),
        m2_factor_mag_sq ~ m2_factor_real^2 + m2_factor_imag^2,

        # Eq. 15.19 - Rayleigh scattering efficiency
        Q_scat ~ (8/3) * α^4 * m2_factor_mag_sq,

        # Eq. 15.21 - Rayleigh absorption efficiency (simplified form for small α)
        Q_abs ~ 4 * α * m2_factor_imag,

        # Eq. 15.4 - Extinction efficiency
        Q_ext ~ Q_scat + Q_abs,

        # Eq. 15.5 - Single-scattering albedo
        ω ~ Q_scat / Q_ext,

        # Mass efficiencies (Eq. 15.41-15.43)
        E_ext ~ 3 * Q_ext / (2 * ρ_p * D_p),
        E_scat ~ 3 * Q_scat / (2 * ρ_p * D_p),
        E_abs ~ 3 * Q_abs / (2 * ρ_p * D_p),
    ]

    return System(eqs, t; name, checks=false)
end

"""
    AerosolExtinction(; name=:AerosolExtinction)

Extinction coefficient calculation for an aerosol population.

Computes the extinction, scattering, and absorption coefficients (b_ext, b_scat, b_abs)
for an ensemble of particles. For a monodisperse population, these are related to the
single-particle cross sections by multiplication with the number concentration.

# Reference
Seinfeld & Pandis (2006), Chapter 15, Equations 15.24-15.30, 15.37-15.40.
"""
@component function AerosolExtinction(; name=:AerosolExtinction)
    @constants begin
        π_const = Float64(π), [description = "Mathematical constant pi (dimensionless)"]
        length_scale = 1.0, [unit = u"m", description = "Length scale for dimensionless conversion"]
    end

    @parameters begin
        D_p = 0.5e-6, [unit = u"m", description = "Particle diameter"]
        λ = 550e-9, [unit = u"m", description = "Wavelength of incident radiation"]
        n_refr = 1.5, [description = "Real part of complex refractive index (dimensionless)"]
        k_refr = 0.0, [description = "Imaginary part of complex refractive index (dimensionless)"]
        ρ_p = 1500.0, [unit = u"kg/m^3", description = "Particle density"]
        N = 1e9, [unit = u"m^-3", description = "Particle number concentration"]
    end

    @variables begin
        α(t), [description = "Size parameter (dimensionless)"]
        Q_ext(t), [description = "Extinction efficiency (dimensionless)"]
        Q_scat(t), [description = "Scattering efficiency (dimensionless)"]
        Q_abs(t), [description = "Absorption efficiency (dimensionless)"]
        b_ext(t), [unit = u"m^-1", description = "Extinction coefficient"]
        b_scat(t), [unit = u"m^-1", description = "Scattering coefficient"]
        b_abs(t), [unit = u"m^-1", description = "Absorption coefficient"]
        ω(t), [description = "Single-scattering albedo (dimensionless)"]
    end

    eqs = [
        # Eq. 15.6 - Size parameter (dimensionless ratio)
        α ~ π_const * (D_p / length_scale) / (λ / length_scale),

        # Mie efficiencies
        Q_ext ~ mie_Q_ext(α, n_refr, k_refr),
        Q_scat ~ mie_Q_scat(α, n_refr, k_refr),
        Q_abs ~ Q_ext - Q_scat,

        # Eq. 15.27 - Extinction coefficient (monodisperse)
        b_ext ~ (π_const / 4) * (D_p / length_scale)^2 * N * Q_ext * length_scale^2,

        # Eq. 15.28 - Scattering and absorption coefficients
        b_scat ~ (π_const / 4) * (D_p / length_scale)^2 * N * Q_scat * length_scale^2,
        b_abs ~ (π_const / 4) * (D_p / length_scale)^2 * N * Q_abs * length_scale^2,

        # Eq. 15.5 - Single-scattering albedo
        ω ~ b_scat / b_ext,
    ]

    return System(eqs, t; name, checks=false)
end

"""
    Visibility(; name=:Visibility)

Visibility calculation from extinction coefficient.

Computes visual range using the Koschmeider equation (Eq. 15.36) and related
visibility metrics. Visual range is the distance at which a black object
against a white background reaches the threshold contrast of 0.02.

# Reference
Seinfeld & Pandis (2006), Chapter 15, Equations 15.31-15.36.
"""
@component function Visibility(; name=:Visibility)
    @constants begin
        # Eq. 15.36 - Koschmeider constant = -ln(0.02) ≈ 3.912
        Koschmeider_const = 3.912, [description = "Koschmeider constant = -ln(0.02) (dimensionless)"]
        # Sea-level Rayleigh scattering coefficient at λ = 520 nm
        b_Rayleigh_sea_level = 13.2e-6, [unit = u"m^-1", description = "Sea-level Rayleigh scattering coefficient at λ=520nm"]
    end

    @parameters begin
        b_ext = 1e-4, [unit = u"m^-1", description = "Total extinction coefficient"]
        b_sg = 13.2e-6, [unit = u"m^-1", description = "Rayleigh scattering coefficient (gas)"]
        b_ag = 0.0, [unit = u"m^-1", description = "Absorption coefficient due to gases (mainly NO₂)"]
    end

    @variables begin
        x_v(t), [unit = u"m", description = "Visual range"]
        τ_per_km(t), [unit = u"km^-1", description = "Optical depth per kilometer"]
        b_sp(t), [unit = u"m^-1", description = "Scattering coefficient due to particles"]
    end

    eqs = [
        # Eq. 15.36 - Koschmeider equation for visual range
        x_v ~ Koschmeider_const / b_ext,

        # Optical depth per km (for reference)
        τ_per_km ~ b_ext * 1000,

        # Eq. 15.30 - Particle scattering coefficient
        # b_scat = b_sg + b_sp, so b_sp = b_ext - b_sg - b_ag (assuming absorption from particles is small)
        b_sp ~ b_ext - b_sg - b_ag,
    ]

    return System(eqs, t; name, checks=false)
end


# =============================================================================
# Mie theory calculation functions (registered for symbolic use)
# =============================================================================
#
# Implementation based on the algorithm from Bohren & Huffman (1983)
# "Absorption and Scattering of Light by Small Particles", Appendix A,
# as referenced in Seinfeld & Pandis (2006) Appendix 15.
#
# Key features for numerical stability:
# - Downward recurrence for logarithmic derivatives A_n
# - Upward recurrence for Riccati-Bessel functions ψ_n and χ_n
# - Careful handling of edge cases

"""
    mie_Q_ext(α, n, k)

Calculate extinction efficiency Q_ext using Mie theory.

Uses the algorithm from Appendix 15 of Seinfeld & Pandis (2006),
based on Bohren & Huffman (1983).

# Arguments
- `α`: Size parameter (= πDₚ/λ)
- `n`: Real part of refractive index
- `k`: Imaginary part of refractive index

# Returns
- Q_ext: Extinction efficiency (dimensionless)
"""
function mie_Q_ext(α::Real, n::Real, k::Real)
    # Handle edge cases
    if α <= 0
        return 0.0
    end
    if α < 0.01
        # Use Rayleigh approximation for very small particles
        return _rayleigh_Q_ext(α, n, k)
    end

    # Compute both Q_ext and Q_scat together for consistency
    Q_ext, _ = _mie_coefficients(α, n, k)
    return Q_ext
end

"""
    mie_Q_scat(α, n, k)

Calculate scattering efficiency Q_scat using Mie theory.

# Arguments
- `α`: Size parameter (= πDₚ/λ)
- `n`: Real part of refractive index
- `k`: Imaginary part of refractive index

# Returns
- Q_scat: Scattering efficiency (dimensionless)
"""
function mie_Q_scat(α::Real, n::Real, k::Real)
    # Handle edge cases
    if α <= 0
        return 0.0
    end
    if α < 0.01
        # Use Rayleigh approximation for very small particles
        return _rayleigh_Q_scat(α, n, k)
    end

    # Compute both Q_ext and Q_scat together for consistency
    _, Q_scat = _mie_coefficients(α, n, k)
    return Q_scat
end

"""
    _mie_coefficients(α, n, k)

Compute Mie extinction and scattering efficiencies together.

Based on the BHMIE algorithm from Bohren & Huffman (1983) "Absorption and
Scattering of Light by Small Particles", Appendix A. This is the standard
reference implementation used widely in atmospheric science.

Returns (Q_ext, Q_scat).
"""
function _mie_coefficients(α::Real, n::Real, k::Real)
    # Complex refractive index (relative to medium, usually air ≈ 1)
    m = Complex(n, k)
    x = α  # Size parameter

    # Number of terms needed (Wiscombe, 1980 criterion)
    nstop = ceil(Int, x + 4.05 * x^(1/3) + 2)

    # For computing logarithmic derivatives of ψ
    nmx = max(nstop + 15, ceil(Int, abs(x * m)) + 15)

    # Compute d_n = d/d(mx)[ln ψ_n(mx)] using downward recurrence
    # This is the logarithmic derivative of ψ_n evaluated at mx
    d = zeros(Complex{Float64}, nmx)
    mx = m * x

    # Starting value (asymptotic approximation)
    for i in nmx:-1:2
        en = i
        d[i-1] = en / mx - 1.0 / (d[i] + en / mx)
    end

    # Initialize pi and tau (angle functions - not needed for Q_ext, Q_scat)
    # Initialize Riccati-Bessel functions at real argument x

    # ψ_0(x) = sin(x), ψ_{-1}(x) = cos(x)
    # χ_0(x) = cos(x), χ_{-1}(x) = -sin(x)
    # Note: ξ_n = ψ_n + i·χ_n

    psi0 = cos(x)
    psi1 = sin(x)
    chi0 = -sin(x)
    chi1 = cos(x)

    Q_ext_sum = 0.0
    Q_scat_sum = 0.0

    for i in 1:nstop
        en = Float64(i)
        fn = (2.0 * en + 1.0) / (en * (en + 1.0))

        # Compute ψ_n and χ_n for order n using upward recurrence
        psi = (2.0 * en - 1.0) * psi1 / x - psi0
        chi = (2.0 * en - 1.0) * chi1 / x - chi0

        xi = Complex(psi, -chi)

        # Compute a_n and b_n (Bohren & Huffman Eq. 4.53)
        # a_n = [ d_n(mx)/m + n/x ] ψ_n(x) - ψ_{n-1}(x)
        #       ----------------------------------------
        #       [ d_n(mx)/m + n/x ] ξ_n(x) - ξ_{n-1}(x)
        #
        # b_n = [ m·d_n(mx) + n/x ] ψ_n(x) - ψ_{n-1}(x)
        #       ----------------------------------------
        #       [ m·d_n(mx) + n/x ] ξ_n(x) - ξ_{n-1}(x)

        dn = d[i]
        xi1 = Complex(psi1, -chi1)

        an_numer = (dn / m + en / x) * psi - psi1
        an_denom = (dn / m + en / x) * xi - xi1
        an = an_numer / an_denom

        bn_numer = (m * dn + en / x) * psi - psi1
        bn_denom = (m * dn + en / x) * xi - xi1
        bn = bn_numer / bn_denom

        # Sum for Q_ext (Eq. 4.62) and Q_scat (Eq. 4.61)
        Q_ext_sum += (2.0 * en + 1.0) * real(an + bn)
        Q_scat_sum += (2.0 * en + 1.0) * (abs2(an) + abs2(bn))

        # Update for next iteration
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
    end

    Q_ext = (2.0 / x^2) * Q_ext_sum
    Q_scat = (2.0 / x^2) * Q_scat_sum

    return Q_ext, Q_scat
end

"""
    _rayleigh_Q_ext(α, n, k)

Rayleigh approximation for extinction efficiency (Eq. 15.18, 15.21).

Valid for α << 1 (particles much smaller than wavelength).
"""
function _rayleigh_Q_ext(α::Real, n::Real, k::Real)
    # (m² - 1)/(m² + 2) where m = n + ik
    m2_minus_1_real, m2_minus_1_imag, m2_plus_2_real, m2_plus_2_imag = _m2_components(n, k)
    denom = m2_plus_2_real^2 + m2_plus_2_imag^2
    imag_factor = (m2_minus_1_imag * m2_plus_2_real - m2_minus_1_real * m2_plus_2_imag) / denom
    mag_sq = (m2_minus_1_real^2 + m2_minus_1_imag^2) / denom

    # Q_ext = Q_abs + Q_scat
    # Q_abs ≈ 4α × Im{(m²-1)/(m²+2)} (Eq. 15.21)
    # Q_scat ≈ (8/3)α⁴ × |(m²-1)/(m²+2)|² (Eq. 15.19)
    return 4 * α * imag_factor + (8/3) * α^4 * mag_sq
end

"""
    _rayleigh_Q_scat(α, n, k)

Rayleigh approximation for scattering efficiency (Eq. 15.19).

Valid for α << 1 (particles much smaller than wavelength).
"""
function _rayleigh_Q_scat(α::Real, n::Real, k::Real)
    m2_minus_1_real, m2_minus_1_imag, m2_plus_2_real, m2_plus_2_imag = _m2_components(n, k)
    denom = m2_plus_2_real^2 + m2_plus_2_imag^2
    mag_sq = (m2_minus_1_real^2 + m2_minus_1_imag^2) / denom

    return (8/3) * α^4 * mag_sq
end

"""
    _m2_components(n, k)

Compute components for (m² - 1) and (m² + 2) where m = n + ik.
Returns (real(m²-1), imag(m²-1), real(m²+2), imag(m²+2)).
"""
function _m2_components(n::Real, k::Real)
    # m² = (n + ik)² = n² - k² + 2ink
    m2_real = n^2 - k^2
    m2_imag = 2 * n * k

    m2_minus_1_real = m2_real - 1
    m2_minus_1_imag = m2_imag

    m2_plus_2_real = m2_real + 2
    m2_plus_2_imag = m2_imag

    return m2_minus_1_real, m2_minus_1_imag, m2_plus_2_real, m2_plus_2_imag
end

"""
    rayleigh_m2_factor_real(n, k)

Real part of (m² - 1)/(m² + 2) for Rayleigh scattering.
"""
function rayleigh_m2_factor_real(n::Real, k::Real)
    m2_minus_1_real, m2_minus_1_imag, m2_plus_2_real, m2_plus_2_imag = _m2_components(n, k)
    denom = m2_plus_2_real^2 + m2_plus_2_imag^2
    return (m2_minus_1_real * m2_plus_2_real + m2_minus_1_imag * m2_plus_2_imag) / denom
end

"""
    rayleigh_m2_factor_imag(n, k)

Imaginary part of (m² - 1)/(m² + 2) for Rayleigh scattering.
"""
function rayleigh_m2_factor_imag(n::Real, k::Real)
    m2_minus_1_real, m2_minus_1_imag, m2_plus_2_real, m2_plus_2_imag = _m2_components(n, k)
    denom = m2_plus_2_real^2 + m2_plus_2_imag^2
    return (m2_minus_1_imag * m2_plus_2_real - m2_minus_1_real * m2_plus_2_imag) / denom
end

# Register the Mie functions for symbolic use in ModelingToolkit
@register_symbolic mie_Q_ext(α::Real, n::Real, k::Real)
@register_symbolic mie_Q_scat(α::Real, n::Real, k::Real)
@register_symbolic rayleigh_m2_factor_real(n::Real, k::Real)
@register_symbolic rayleigh_m2_factor_imag(n::Real, k::Real)
