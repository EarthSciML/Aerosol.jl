# Aerosol Size Distribution

## Overview

This module implements the multi-modal lognormal aerosol size distribution
parameterization described in Chapter 8 ("Properties of the Atmospheric Aerosol")
of Seinfeld and Pandis (2006).

The aerosol particle size distribution is represented as a sum of lognormal modes
(Eq. 8.54), where each mode is characterized by a number concentration ``N_i``,
geometric median diameter ``D_{g,i}``, and geometric standard deviation ``\sigma_{g,i}``.
This parameterization is widely used to represent the number, surface area, and volume
distributions of atmospheric aerosol populations across different environments.

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and
Physics: From Air Pollution to Climate Change*, 2nd Edition, John Wiley & Sons,
Chapter 8, pp. 350–395.

```@docs
AerosolDistribution
UrbanAerosol
MarineAerosol
RuralAerosol
RemoteContinentalAerosol
FreeTroposphereAerosol
PolarAerosol
DesertAerosol
```

## Implementation

The `AerosolDistribution` component represents a multi-modal lognormal distribution
with the following key equations:

- **Number distribution** (Eq. 8.54): ``\frac{dN}{d\log D_p} = \sum_i \frac{N_i}{\sqrt{2\pi}\,\log\sigma_{g,i}} \exp\!\left(-\frac{(\log D_p - \log D_{g,i})^2}{2\log^2\sigma_{g,i}}\right)``
- **Surface area distribution** (Eq. 8.50): Same form with surface area median diameter ``D_{s,i} = D_{g,i}\exp(2\ln^2\sigma_{g,i})``
- **Volume distribution** (Eq. 8.51): Same form with volume median diameter ``D_{v,i} = D_{g,i}\exp(3\ln^2\sigma_{g,i})``
- **Vertical mass profile** (Eq. 8.55): ``M(z) = M(0)\exp(-z/H_p)``

Seven predefined aerosol distribution types are provided from Table 8.3, each with
three lognormal modes calibrated to represent typical atmospheric conditions.

### State Variables

```@example size_dist
using DataFrames, ModelingToolkit, Aerosol, DynamicQuantities

sys = AerosolDistribution(3)
vars = unknowns(sys)
DataFrame(
    :Name => [string(v) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example size_dist
pars = parameters(sys)
DataFrame(
    :Name => [string(p) for p in pars],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in pars],
    :Description => [ModelingToolkit.getdescription(p) for p in pars]
)
```

### Equations

```@example size_dist
eqs = equations(sys)
```

## Analysis

### Urban Aerosol Distribution (Figure 8.14)

The urban aerosol distribution from Table 8.3 shows three distinct modes:
a nucleation mode, an Aitken mode, and an accumulation mode.

```@example size_dist
using Plots

urban = UrbanAerosol()
defs = ModelingToolkit.defaults(urban)

# Extract mode parameters from defaults
N_modes = Float64[]
D_g_modes = Float64[]
logσ_modes = Float64[]
for i in 1:3
    n_key = filter(k -> string(k) == "UrbanAerosol₊N[$i]", collect(keys(defs)))
    d_key = filter(k -> string(k) == "UrbanAerosol₊D_g[$i]", collect(keys(defs)))
    s_key = filter(k -> string(k) == "UrbanAerosol₊logσ[$i]", collect(keys(defs)))
    push!(N_modes, defs[n_key[1]])
    push!(D_g_modes, defs[d_key[1]])
    push!(logσ_modes, defs[s_key[1]])
end

# Compute distribution over a range of diameters
D_p_range = 10 .^ range(log10(1e-9), log10(100e-6), length=500)

function lognormal_mode(D_p, N, D_g, logσ)
    return N / (sqrt(2π) * logσ) * exp(-(log10(D_p / D_g))^2 / (2 * logσ^2))
end

n_N_total = zeros(length(D_p_range))
mode_contributions = [zeros(length(D_p_range)) for _ in 1:3]

for (j, D_p) in enumerate(D_p_range)
    for i in 1:3
        val = lognormal_mode(D_p, N_modes[i], D_g_modes[i], logσ_modes[i])
        mode_contributions[i][j] = val
        n_N_total[j] += val
    end
end

# Convert to cm^-3 for plotting (standard convention)
p = plot(D_p_range * 1e6, n_N_total / 1e6,
    xscale=:log10, xlabel="Particle Diameter (μm)",
    ylabel="dN/d(log Dp) (cm⁻³)",
    label="Total", linewidth=2, color=:black,
    title="Urban Aerosol Size Distribution (Table 8.3)")
for i in 1:3
    plot!(p, D_p_range * 1e6, mode_contributions[i] / 1e6,
        label="Mode $i", linestyle=:dash)
end
p
```

### Comparison of All Aerosol Types (Figures 8.14–8.20)

The following figure compares the number distributions for all seven
aerosol environment types from Table 8.3.

```@example size_dist
# Define all distribution types and their parameters from Table 8.3
dist_params = Dict(
    "Urban" => [(9.93e4, 0.013e-6, 0.245), (1.11e3, 0.014e-6, 0.666), (3.64e4, 0.050e-6, 0.337)],
    "Marine" => [(133.0, 0.008e-6, 0.657), (66.6, 0.266e-6, 0.210), (3.06, 0.580e-6, 0.396)],
    "Rural" => [(6.65e3, 0.015e-6, 0.225), (147.0, 0.054e-6, 0.557), (1990.0, 0.084e-6, 0.266)],
    "Remote continental" => [(3200.0, 0.020e-6, 0.161), (2900.0, 0.116e-6, 0.217), (0.300, 1.800e-6, 0.380)],
    "Free troposphere" => [(129.0, 0.007e-6, 0.645), (59.7, 0.250e-6, 0.253), (63.5, 0.520e-6, 0.425)],
    "Polar" => [(21.7, 0.138e-6, 0.164), (0.186, 0.750e-6, 0.521), (3.04e-4, 8.600e-6, 0.420)],
    "Desert" => [(726.0, 0.002e-6, 0.247), (114.0, 0.038e-6, 0.770), (0.178, 21.60e-6, 0.438)],
)

D_p_range = 10 .^ range(log10(1e-9), log10(100e-6), length=500)
p = plot(xlabel="Particle Diameter (μm)", ylabel="dN/d(log Dp) (cm⁻³)",
    xscale=:log10, yscale=:log10, legend=:topright,
    title="Model Aerosol Size Distributions (Table 8.3)",
    ylims=(1e-3, 1e6))

for (env_name, modes) in sort(collect(dist_params), by=x->x[1])
    n_N = zeros(length(D_p_range))
    for (j, D_p) in enumerate(D_p_range)
        for (N, D_g, logσ) in modes
            # N is in cm^-3 here (not SI), compute directly in cm^-3
            n_N[j] += lognormal_mode(D_p, N, D_g, logσ)
        end
    end
    plot!(p, D_p_range * 1e6, n_N, label=env_name, linewidth=1.5)
end
p
```

### Volume Distributions

The volume distribution emphasizes larger particles and is important for
understanding aerosol mass loading. For each mode, the volume distribution
peaks at the volume median diameter ``D_{v,i} = D_{g,i}\exp(3\ln^2\sigma_{g,i})``
(Eq. 8.52).

```@example size_dist
# Urban volume distribution
p = plot(xlabel="Particle Diameter (μm)", ylabel="dV/d(log Dp) (μm³/cm³)",
    xscale=:log10, title="Urban Aerosol Volume Distribution",
    legend=:topright)

# Volume distribution: n_V^o = (π/6) * N * D_pv^3 / (sqrt(2π) * logσ) * exp(...)
# where D_pv = D_g * exp(3 * (logσ * ln10)^2)
n_V_total = zeros(length(D_p_range))
for (j, D_p) in enumerate(D_p_range)
    for i in 1:3
        N = [9.93e4, 1.11e3, 3.64e4][i]  # cm^-3
        D_g = [0.013e-6, 0.014e-6, 0.050e-6][i]
        logσ = [0.245, 0.666, 0.337][i]
        lnσ = logσ * log(10)
        D_pv = D_g * exp(3 * lnσ^2)
        # Eq. 8.51 - volume distribution in log-space
        n_V = (π / 6) * N * D_pv^3 / (sqrt(2π) * logσ) *
              exp(-(log10(D_p / D_pv))^2 / (2 * logσ^2))
        n_V_total[j] += n_V
    end
end

# Convert from m³/cm³ to μm³/cm³ (multiply by 1e18)
plot!(p, D_p_range * 1e6, n_V_total * 1e18, label="Total", linewidth=2, color=:black)
p
```

### Vertical Mass Profile (Eq. 8.55)

The vertical variation of aerosol mass concentration follows an exponential
decay with a characteristic scale height ``H_p``.

```@example size_dist
z_range = range(0, 10000, length=100)
H_p_values = [500, 1000, 2000, 5000]

p = plot(xlabel="M(z)/M(0)", ylabel="Altitude (m)",
    title="Vertical Aerosol Mass Profile (Eq. 8.55)",
    legend=:topright)

for H_p in H_p_values
    M_ratio = exp.(-z_range / H_p)
    plot!(p, M_ratio, z_range, label="Hp = $H_p m", linewidth=1.5)
end
p
```
