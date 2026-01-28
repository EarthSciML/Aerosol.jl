# Interaction of Aerosols with Radiation

## Overview

This module implements the theory of light scattering, absorption, and extinction by
spherical aerosol particles, based on Mie theory and its limiting cases. These optical
properties are fundamental to understanding visibility, radiative transfer, and the
direct radiative forcing of aerosols.

**Reference**: Seinfeld, J. H., & Pandis, S. N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change* (2nd ed.), Chapter 15: Interaction of Aerosols with
Radiation.

```@docs
MieScattering
```

```@docs
RayleighScattering
```

```@docs
AerosolExtinction
```

```@docs
Visibility
```

## Implementation

### MieScattering Component

The `MieScattering` component computes single-particle optical properties using full Mie
theory. The key equations are:

- **Size parameter** (Eq. 15.6): ``\alpha = \pi D_p / \lambda``
- **Extinction efficiency** (Eq. 15.14): ``Q_{\text{ext}} = \frac{2}{\alpha^2} \sum_{k=1}^{\infty} (2k+1) \text{Re}(a_k + b_k)``
- **Scattering efficiency** (Eq. 15.13): ``Q_{\text{scat}} = \frac{2}{\alpha^2} \sum_{k=1}^{\infty} (2k+1)(|a_k|^2 + |b_k|^2)``
- **Absorption efficiency** (Eq. 15.4): ``Q_{\text{abs}} = Q_{\text{ext}} - Q_{\text{scat}}``
- **Single-scattering albedo** (Eq. 15.5): ``\omega = Q_{\text{scat}} / Q_{\text{ext}}``
- **Mass extinction efficiency** (Eq. 15.41): ``E_{\text{ext}} = \frac{3 Q_{\text{ext}}}{2 \rho_p D_p}``

where ``a_k`` and ``b_k`` are the Mie coefficients computed using the Bohren & Huffman (1983)
algorithm.

#### State Variables

```@example mie
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Aerosol

sys = MieScattering()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Parameters

```@example mie
ps = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in ps],
    :Description => [ModelingToolkit.getdescription(p) for p in ps]
)
```

#### Equations

```@example mie
equations(sys)
```

### RayleighScattering Component

The `RayleighScattering` component provides the closed-form Rayleigh approximation valid
when ``\alpha \ll 1`` (particles much smaller than the wavelength):

- **Rayleigh scattering efficiency** (Eq. 15.19): ``Q_{\text{scat}} = \frac{8}{3}\alpha^4 \left|\frac{m^2-1}{m^2+2}\right|^2``
- **Rayleigh absorption efficiency** (Eq. 15.21): ``Q_{\text{abs}} = 4\alpha \, \text{Im}\left\{\frac{m^2-1}{m^2+2}\right\}``

### AerosolExtinction Component

Computes the extinction, scattering, and absorption coefficients for an aerosol population:

- **Extinction coefficient** (Eq. 15.27): ``b_{\text{ext}} = \frac{\pi D_p^2}{4} N Q_{\text{ext}}``

### Visibility Component

Uses the Koschmeider equation to compute visual range:

- **Visual range** (Eq. 15.36): ``x_v = \frac{3.912}{b_{\text{ext}}}``

## Analysis

### Mie Efficiency vs. Size Parameter (Fig. 15.4)

The extinction efficiency ``Q_{\text{ext}}`` oscillates as a function of the size parameter
``\alpha``, converging to 2 for large particles (the "extinction paradox", Eq. 15.23).

```@example mie
using Plots

alpha_range = 0.01:0.05:30.0
n = 1.5  # Real part of refractive index
k = 0.0  # Non-absorbing

Q_ext_vals = [Aerosol.mie_Q_ext(α, n, k) for α in alpha_range]
Q_scat_vals = [Aerosol.mie_Q_scat(α, n, k) for α in alpha_range]

plot(alpha_range, Q_ext_vals, label="Q_ext", xlabel="Size parameter α",
     ylabel="Efficiency", title="Mie Extinction Efficiency (n=$n, k=$k)",
     linewidth=2, legend=:topright)
hline!([2.0], linestyle=:dash, label="Q_ext → 2 (extinction paradox)", color=:gray)
```

### Effect of Absorption on Scattering (Fig. 15.5)

For absorbing particles (``k > 0``), the single-scattering albedo ``\omega`` decreases,
indicating a larger fraction of extinction is due to absorption.

```@example mie
alpha_range = 0.01:0.05:20.0

# Non-absorbing (water-like)
Q_ext_water = [Aerosol.mie_Q_ext(α, 1.33, 0.0) for α in alpha_range]
Q_scat_water = [Aerosol.mie_Q_scat(α, 1.33, 0.0) for α in alpha_range]
Q_abs_water = Q_ext_water .- Q_scat_water

# Moderately absorbing
Q_ext_mod = [Aerosol.mie_Q_ext(α, 1.5, 0.1) for α in alpha_range]
Q_scat_mod = [Aerosol.mie_Q_scat(α, 1.5, 0.1) for α in alpha_range]
Q_abs_mod = Q_ext_mod .- Q_scat_mod

# Strongly absorbing (carbon-like)
Q_ext_carbon = [Aerosol.mie_Q_ext(α, 1.95, 0.79) for α in alpha_range]
Q_scat_carbon = [Aerosol.mie_Q_scat(α, 1.95, 0.79) for α in alpha_range]
Q_abs_carbon = Q_ext_carbon .- Q_scat_carbon

p1 = plot(alpha_range, Q_ext_water, label="Q_ext (water)", linewidth=2)
plot!(alpha_range, Q_ext_mod, label="Q_ext (n=1.5, k=0.1)", linewidth=2)
plot!(alpha_range, Q_ext_carbon, label="Q_ext (carbon)", linewidth=2)
xlabel!("Size parameter α")
ylabel!("Q_ext")
title!("Extinction Efficiency for Different Materials")

p2 = plot(alpha_range, Q_abs_water, label="Q_abs (water)", linewidth=2)
plot!(alpha_range, Q_abs_mod, label="Q_abs (n=1.5, k=0.1)", linewidth=2)
plot!(alpha_range, Q_abs_carbon, label="Q_abs (carbon)", linewidth=2)
xlabel!("Size parameter α")
ylabel!("Q_abs")
title!("Absorption Efficiency for Different Materials")

plot(p1, p2, layout=(2,1), size=(700, 600))
```

### Single-Scattering Albedo vs. Size Parameter

The single-scattering albedo ``\omega = Q_{\text{scat}} / Q_{\text{ext}}`` characterizes
the relative importance of scattering vs. absorption.

```@example mie
alpha_range = 0.1:0.05:20.0

ω_water = [Aerosol.mie_Q_scat(α, 1.33, 0.0) / Aerosol.mie_Q_ext(α, 1.33, 0.0) for α in alpha_range]
ω_mod = [Aerosol.mie_Q_scat(α, 1.5, 0.1) / Aerosol.mie_Q_ext(α, 1.5, 0.1) for α in alpha_range]
ω_carbon = [Aerosol.mie_Q_scat(α, 1.95, 0.79) / Aerosol.mie_Q_ext(α, 1.95, 0.79) for α in alpha_range]

plot(alpha_range, ω_water, label="Water (n=1.33, k≈0)", linewidth=2)
plot!(alpha_range, ω_mod, label="Moderate absorption (n=1.5, k=0.1)", linewidth=2)
plot!(alpha_range, ω_carbon, label="Carbon (n=1.95, k=0.79)", linewidth=2)
xlabel!("Size parameter α")
ylabel!("Single-scattering albedo ω")
title!("Single-Scattering Albedo (Eq. 15.5)")
ylims!(0, 1.05)
```

### Mass Extinction Efficiency (Fig. 15.7)

The mass extinction efficiency ``E_{\text{ext}}`` peaks at particle diameters comparable
to the wavelength of light, making submicron particles most efficient at scattering
and absorbing radiation per unit mass.

```@example mie
D_p_range = 0.01e-6:0.01e-6:2.5e-6  # 0.01 to 2.5 μm
λ = 550e-9  # 550 nm

# Ammonium sulfate: n=1.521, k≈0, ρ=1770 kg/m³
n_sulfate = 1.521
k_sulfate = 0.0
ρ_sulfate = 1770.0

E_ext_sulfate = Float64[]
for D_p in D_p_range
    α = π * D_p / λ
    Q_ext = Aerosol.mie_Q_ext(α, n_sulfate, k_sulfate)
    push!(E_ext_sulfate, 3 * Q_ext / (2 * ρ_sulfate * D_p))
end

# Carbon: n=1.95, k=0.79, ρ=2000 kg/m³
n_carbon = 1.95
k_carbon = 0.79
ρ_carbon = 2000.0

E_ext_carbon = Float64[]
for D_p in D_p_range
    α = π * D_p / λ
    Q_ext = Aerosol.mie_Q_ext(α, n_carbon, k_carbon)
    push!(E_ext_carbon, 3 * Q_ext / (2 * ρ_carbon * D_p))
end

D_p_um = collect(D_p_range) .* 1e6  # Convert to μm for plotting

plot(D_p_um, E_ext_sulfate, label="(NH₄)₂SO₄", linewidth=2)
plot!(D_p_um, E_ext_carbon, label="Elemental carbon", linewidth=2)
xlabel!("Particle diameter Dₚ (μm)")
ylabel!("Mass extinction efficiency E_ext (m²/kg)")
title!("Mass Extinction Efficiency at λ=550 nm (Eq. 15.41)")
```

### Rayleigh Scattering Regime

In the Rayleigh regime (``\alpha \ll 1``), the full Mie solution reduces to simple
closed-form expressions. The key feature is that scattering efficiency scales as
``\alpha^4`` (equivalently ``1/\lambda^4``), explaining why the sky is blue.

```@example mie
alpha_range_small = 0.001:0.001:0.1

Q_scat_mie = [Aerosol.mie_Q_scat(α, 1.5, 0.0) for α in alpha_range_small]
Q_scat_rayleigh = [(8/3) * α^4 * ((1.5^2 - 1)/(1.5^2 + 2))^2 for α in alpha_range_small]

plot(alpha_range_small, Q_scat_mie, label="Mie (full)", linewidth=2)
plot!(alpha_range_small, Q_scat_rayleigh, label="Rayleigh (Eq. 15.19)", linewidth=2, linestyle=:dash)
xlabel!("Size parameter α")
ylabel!("Q_scat")
title!("Rayleigh vs. Mie Scattering (n=1.5, k=0)")
```

### Visibility and Extinction

The Koschmeider equation (Eq. 15.36) relates visual range to extinction coefficient.
At sea level in a clean Rayleigh atmosphere (``b_{\text{ext}} \approx 13.2 \times 10^{-6} \text{ m}^{-1}``),
the maximum visual range is approximately 296 km.

```@example mie
b_ext_range = 1e-5:1e-5:1e-3  # Extinction coefficient range (m^-1)
x_v = 3.912 ./ b_ext_range  # Koschmeider equation

plot(b_ext_range .* 1e6, x_v ./ 1000, label="Visual range",
     linewidth=2, xscale=:log10)
xlabel!("Extinction coefficient b_ext (10⁻⁶ m⁻¹)")
ylabel!("Visual range xᵥ (km)")
title!("Koschmeider Equation (Eq. 15.36)")
```
