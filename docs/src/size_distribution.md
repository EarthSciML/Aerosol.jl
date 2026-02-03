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
using DataFrames, ModelingToolkit, Symbolics, Aerosol, DynamicQuantities

sys = AerosolDistribution(3)
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example size_dist
pars = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in pars],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in pars],
    :Description => [ModelingToolkit.getdescription(p) for p in pars]
)
```

### Equations

```@example size_dist
eqs = equations(sys)
```

## Analysis

```@example size_dist
# Helper to get the base name of a symbolic variable (part after last ₊),
# and strip the (t) suffix if present.
function get_base_name(sym)
    name = string(Symbolics.tosymbol(sym, escape = false))
    name = contains(name, "₊") ? split(name, "₊")[end] : name
    name = replace(name, r"\(t\)$" => "")
    return name
end

# Helper to evaluate a system equation's RHS by substituting parameter values.
# Handles namespaced variables by matching the actual symbols in the equation RHS.
function eval_eq(sys, var, param_vals)
    var_base = get_base_name(var)

    # Find the equation for this variable
    matching_eqs = filter(equations(sys)) do eq
        lhs_base = get_base_name(eq.lhs)
        return lhs_base == var_base
    end
    eq = only(matching_eqs)

    # Get all variables that appear in the RHS
    rhs_vars = Symbolics.get_variables(eq.rhs)

    # Build substitution dict by matching RHS variables with:
    # 1. System defaults (for constants like π_c, ln10)
    # 2. User-provided values (matching by string representation)
    subs = Dict{Any, Any}()
    defs = ModelingToolkit.defaults(sys)

    for rv in rhs_vars
        rv_str = string(rv)

        # Check system defaults
        for (k, v) in defs
            k_str = string(k)
            if rv_str == k_str || endswith(k_str, "₊" * rv_str)
                subs[rv] = v
            end
        end

        # Check user-provided values (these override defaults)
        for (k, v) in param_vals
            k_str = string(k)
            if rv_str == k_str || endswith(k_str, "₊" * rv_str)
                subs[rv] = v
            end
        end
    end

    result = Symbolics.substitute(eq.rhs, subs)
    return Float64(result)
end
nothing # hide
```

### Urban Aerosol Distribution (Figure 8.14)

The urban aerosol distribution from Table 8.3 shows three distinct modes:
a nucleation mode, an Aitken mode, and an accumulation mode.

```@example size_dist
using Plots

urban = UrbanAerosol()

# Evaluate the system's n_N_o equation over a range of diameters
D_p_range = 10 .^ range(log10(1e-9), log10(100e-6), length = 500)
n_N_total = map(D_p_range) do D_p_val
    eval_eq(urban, urban.n_N_o, Dict(urban.D_p => D_p_val))
end

# Also evaluate individual mode contributions using single-mode systems
mode_contributions = [zeros(length(D_p_range)) for _ in 1:3]
defs = ModelingToolkit.defaults(urban)
for i in 1:3
    n_key = only(filter(k -> contains(string(k), "N[$i]"), collect(keys(defs))))
    d_key = only(filter(k -> contains(string(k), "D_g[$i]"), collect(keys(defs))))
    s_key = only(filter(k -> contains(string(k), "logσ[$i]"), collect(keys(defs))))
    mode_sys = AerosolDistribution(1; name = Symbol("mode_$i"))
    for (j, D_p_val) in enumerate(D_p_range)
        mode_contributions[i][j] = eval_eq(mode_sys,
            mode_sys.n_N_o,
            Dict(
                mode_sys.N[1] => defs[n_key],
                mode_sys.D_g[1] => defs[d_key],
                mode_sys.logσ[1] => defs[s_key],
                mode_sys.D_p => D_p_val
            ))
    end
end

# Convert to cm^-3 for plotting (standard convention)
p = plot(D_p_range * 1e6, n_N_total / 1e6,
    xscale = :log10, xlabel = "Particle Diameter (μm)",
    ylabel = "dN/d(log Dp) (cm⁻³)",
    label = "Total", linewidth = 2, color = :black,
    title = "Urban Aerosol Size Distribution (Table 8.3)")
for i in 1:3
    plot!(p, D_p_range * 1e6, mode_contributions[i] / 1e6,
        label = "Mode $i", linestyle = :dash)
end
p
```

### Comparison of All Aerosol Types (Figures 8.14–8.20)

The following figure compares the number distributions for all seven
aerosol environment types from Table 8.3.

```@example size_dist
dist_constructors = [
    ("Desert", DesertAerosol),
    ("Free troposphere", FreeTroposphereAerosol),
    ("Marine", MarineAerosol),
    ("Polar", PolarAerosol),
    ("Remote continental", RemoteContinentalAerosol),
    ("Rural", RuralAerosol),
    ("Urban", UrbanAerosol)
]

D_p_range = 10 .^ range(log10(1e-9), log10(100e-6), length = 500)
p = plot(xlabel = "Particle Diameter (μm)", ylabel = "dN/d(log Dp) (cm⁻³)",
    xscale = :log10, yscale = :log10, legend = :topright,
    title = "Model Aerosol Size Distributions (Table 8.3)",
    ylims = (1e-3, 1e6))

for (env_name, constructor) in dist_constructors
    dist_sys = constructor()
    n_N = map(D_p_range) do D_p_val
        eval_eq(dist_sys, dist_sys.n_N_o, Dict(dist_sys.D_p => D_p_val))
    end
    # Convert from m^-3 to cm^-3 for plotting
    plot!(p, D_p_range * 1e6, n_N / 1e6, label = env_name, linewidth = 1.5)
end
p
```

### Volume Distributions

The volume distribution emphasizes larger particles and is important for
understanding aerosol mass loading. For each mode, the volume distribution
peaks at the volume median diameter ``D_{v,i} = D_{g,i}\exp(3\ln^2\sigma_{g,i})``
(Eq. 8.52).

```@example size_dist
# Urban volume distribution using the system's n_V_o equation
urban = UrbanAerosol()
D_p_range = 10 .^ range(log10(1e-9), log10(100e-6), length = 500)

n_V_total = map(D_p_range) do D_p_val
    eval_eq(urban, urban.n_V_o, Dict(urban.D_p => D_p_val))
end

p = plot(xlabel = "Particle Diameter (μm)", ylabel = "dV/d(log Dp) (μm³/cm³)",
    xscale = :log10, title = "Urban Aerosol Volume Distribution",
    legend = :topright)

# Convert from m³/m³ to μm³/cm³: m³→μm³ is ×1e18, m⁻³→cm⁻³ is ÷1e6, net ×1e12
plot!(p, D_p_range * 1e6, n_V_total * 1e12, label = "Total", linewidth = 2, color = :black)
p
```

### Vertical Mass Profile (Eq. 8.55)

The vertical variation of aerosol mass concentration follows an exponential
decay with a characteristic scale height ``H_p``.

```@example size_dist
# Evaluate the system's M_z equation for different scale heights
sys = AerosolDistribution(3)
z_range = range(0, 10000, length = 100)
H_p_values = [500, 1000, 2000, 5000]

base_params = Dict(
    sys.N[1] => 1e11, sys.D_g[1] => 1e-7, sys.logσ[1] => 0.3,
    sys.N[2] => 1e10, sys.D_g[2] => 1e-6, sys.logσ[2] => 0.3,
    sys.N[3] => 1e9, sys.D_g[3] => 1e-5, sys.logσ[3] => 0.3,
    sys.D_p => 1e-7
)

p = plot(xlabel = "M(z)/M(0)", ylabel = "Altitude (m)",
    title = "Vertical Aerosol Mass Profile (Eq. 8.55)",
    legend = :topright)

for H_p in H_p_values
    M_ratio = map(z_range) do z_val
        eval_eq(sys, sys.M_z, merge(base_params, Dict(sys.z => z_val, sys.H_p => Float64(H_p))))
    end
    plot!(p, M_ratio, z_range, label = "Hp = $H_p m", linewidth = 1.5)
end
p
```
