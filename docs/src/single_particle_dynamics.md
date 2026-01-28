# Single Particle Dynamics

## Overview

This module implements equations for the dynamics of single aerosol particles, covering the fundamental physics that govern particle motion, transport, and interactions with gas molecules. These equations are essential for understanding aerosol behavior in the atmosphere, including sedimentation, diffusion, and electrical migration.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 2nd Edition, Chapter 9: Dynamics of Single Aerosol Particles, pp. 396-433. John Wiley & Sons, Inc.

```@docs
SingleParticleDynamics
MeanFreePath
SlipCorrection
SettlingVelocity
BrownianDiffusion
ParticleMobility
ElectricalMobility
StokesNumber
AerodynamicDiameter
```

## Implementation

The implementation provides several components that can be used individually or combined through the comprehensive `SingleParticleDynamics` component.

### Key Physical Quantities

The implementation calculates the following key quantities:

| Quantity | Symbol | Equation | Description |
|----------|--------|----------|-------------|
| Mean free path | λ | Eq. 9.6 | Average distance traveled by air molecules between collisions |
| Knudsen number | Kn | Eq. 9.1 | Ratio of mean free path to particle diameter |
| Slip correction | Cc | Eq. 9.34 | Correction for non-continuum effects |
| Relaxation time | τ | Eq. 9.38 | Time constant for velocity adjustment |
| Terminal velocity | vt | Eq. 9.42 | Steady-state settling velocity |
| Diffusion coefficient | D | Eq. 9.73 | Brownian diffusivity |
| Particle mobility | B | Eq. 9.78 | Mechanical mobility |

### State Variables

```@example spd
using Aerosol, ModelingToolkit, DataFrames, DynamicQuantities

sys = SingleParticleDynamics()
vars = unknowns(sys)
DataFrame(
    :Name => [string(v) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example spd
params = parameters(sys)
DataFrame(
    :Name => [string(p) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example spd
equations(sys)
```

## Analysis

### Slip Correction Factor vs Particle Size

The Cunningham slip correction factor (Eq. 9.34) accounts for non-continuum effects when particle size approaches the mean free path of gas molecules. For particles much larger than the mean free path, Cc approaches 1 (continuum regime). For smaller particles, Cc increases significantly.

```@example spd
using Plots, NonlinearSolve

sys = SingleParticleDynamics()
csys = mtkcompile(sys)

# Calculate Cc for different particle sizes
D_p_range = 10 .^ range(-9, -4, length=50)  # 1 nm to 100 μm
C_c_values = Float64[]
Kn_values = Float64[]

for D_p in D_p_range
    prob = NonlinearProblem(csys, Dict(), Dict(csys.D_p => D_p))
    sol = solve(prob)
    push!(C_c_values, sol[csys.C_c])
    push!(Kn_values, sol[csys.Kn])
end

p1 = plot(D_p_range * 1e6, C_c_values,
    xscale=:log10, yscale=:log10,
    xlabel="Particle Diameter (μm)",
    ylabel="Slip Correction Factor Cc",
    title="Cunningham Slip Correction (Eq. 9.34)",
    legend=false, linewidth=2)
hline!([1.0], linestyle=:dash, color=:gray, label="Continuum limit")
```

### Terminal Settling Velocity vs Particle Size

The terminal settling velocity (Eq. 9.42) shows the balance between gravitational force and drag. Larger particles settle faster due to their greater mass.

```@example spd
v_t_values = Float64[]

for D_p in D_p_range
    prob = NonlinearProblem(csys, Dict(), Dict(csys.D_p => D_p))
    sol = solve(prob)
    push!(v_t_values, sol[csys.v_t])
end

p2 = plot(D_p_range * 1e6, v_t_values * 100 * 3600,  # Convert to cm/h
    xscale=:log10, yscale=:log10,
    xlabel="Particle Diameter (μm)",
    ylabel="Terminal Velocity (cm/h)",
    title="Settling Velocity (Eq. 9.42)",
    legend=false, linewidth=2)
```

### Brownian Diffusion Coefficient vs Particle Size

The Stokes-Einstein-Sutherland relation (Eq. 9.73) shows that smaller particles have larger diffusion coefficients, allowing them to spread more rapidly through random Brownian motion.

```@example spd
D_B_values = Float64[]

for D_p in D_p_range
    prob = NonlinearProblem(csys, Dict(), Dict(csys.D_p => D_p))
    sol = solve(prob)
    push!(D_B_values, sol[csys.D_B])
end

p3 = plot(D_p_range * 1e6, D_B_values * 1e4,  # Convert to cm²/s
    xscale=:log10, yscale=:log10,
    xlabel="Particle Diameter (μm)",
    ylabel="Diffusion Coefficient (cm²/s)",
    title="Brownian Diffusivity (Eq. 9.73)",
    legend=false, linewidth=2)
```

### Comparison of Settling and Diffusion Timescales

The relative importance of gravitational settling versus Brownian diffusion depends on particle size. The characteristic timescale ratio τ_ds = 4D/vt² (Eq. 9.85) indicates which process dominates:

```@example spd
tau_ds_values = 4 .* D_B_values ./ (v_t_values.^2)

p4 = plot(D_p_range * 1e6, tau_ds_values,
    xscale=:log10, yscale=:log10,
    xlabel="Particle Diameter (μm)",
    ylabel="Characteristic Time τ_ds (s)",
    title="Diffusion-Settling Timescale (Eq. 9.85)",
    legend=false, linewidth=2)
hline!([1.0], linestyle=:dash, color=:gray)
```

### Summary Plot

```@example spd
plot(p1, p2, p3, p4, layout=(2,2), size=(800, 600))
```

### Validation Against Table 9.3: Slip Correction Factors

The following compares computed slip correction factors against values from Table 9.3 in Seinfeld & Pandis (2006):

```@example spd
# Reference values from Table 9.3
table_data = [
    (0.001, 216.0),
    (0.01, 22.2),
    (0.1, 2.85),
    (1.0, 1.164),
    (10.0, 1.016),
    (100.0, 1.002)
]

D_p_table = [d[1] * 1e-6 for d in table_data]  # Convert μm to m
Cc_table = [d[2] for d in table_data]

# Computed values
Cc_computed = Float64[]
for D_p in D_p_table
    prob = NonlinearProblem(csys, Dict(), Dict(csys.D_p => D_p))
    sol = solve(prob)
    push!(Cc_computed, sol[csys.C_c])
end

DataFrame(
    Symbol("Dp (μm)") => [d[1] for d in table_data],
    Symbol("Cc (Table 9.3)") => Cc_table,
    Symbol("Cc (Computed)") => round.(Cc_computed, digits=2),
    Symbol("Error (%)") => round.(abs.(Cc_computed .- Cc_table) ./ Cc_table .* 100, digits=1)
)
```

### Flow Regime Classification

The Knudsen number determines the applicable flow regime:

| Regime | Kn Range | Description |
|--------|----------|-------------|
| Continuum | Kn → 0 | Stokes' law applies (Cc ≈ 1) |
| Transition | Kn ~ 1 | Slip correction required |
| Free molecular | Kn → ∞ | Kinetic theory applies |

```@example spd
plot(D_p_range * 1e6, Kn_values,
    xscale=:log10, yscale=:log10,
    xlabel="Particle Diameter (μm)",
    ylabel="Knudsen Number",
    title="Flow Regime Classification",
    legend=false, linewidth=2)
hline!([0.1, 10], linestyle=:dash, color=:gray, label="Regime boundaries")
annotate!([(0.001, 0.05, text("Continuum", 8)),
           (0.1, 5, text("Transition", 8)),
           (10, 500, text("Free Molecular", 8))])
```
