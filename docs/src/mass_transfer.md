# Mass Transfer to Atmospheric Particles

## Overview

This module implements the fundamental equations describing mass transfer between gas and particle phases in the atmosphere. The treatment covers three transport regimes based on the Knudsen number (Kn = λ/Rp, ratio of mean free path to particle radius):

  - **Continuum regime** (Kn ≪ 1): Applicable to large particles where Fick's law of diffusion governs transport
  - **Kinetic regime** (Kn ≫ 1): Applicable to small particles where molecular kinetic theory governs transport
  - **Transition regime** (Kn ~ 1): Requires interpolation formulas bridging the two limiting cases

The implementation includes:

  - Transition regime correction factors (Fuchs-Sutugin, Dahneke)
  - Maxwellian flux for particle growth/evaporation
  - Mass transfer coefficients for gas-particle exchange
  - Characteristic timescales for various mass transport processes
  - Aqueous-phase chemistry and mass transport coupling

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006) "Atmospheric Chemistry and Physics: From Air Pollution to Climate Change", 2nd Edition, John Wiley & Sons, Chapter 12, pp. 537-587.

```@docs
MassTransfer
FuchsSutugin
Dahneke
MaxwellianFlux
MassTransferCoefficient
UptakeCoefficient
```

## Implementation

### Transition Regime Mass Transfer

The Fuchs-Sutugin formula is the most widely used expression for mass transfer across the full range of Knudsen numbers:

```@example mass_transfer
using Aerosol
using ModelingToolkit
import ModelingToolkit: mtkcompile
using NonlinearSolve
using Plots

# Create the Fuchs-Sutugin system
sys = FuchsSutugin()
compiled = mtkcompile(sys)

# Calculate correction factor across Knudsen number range
Kn_range = 10 .^ range(-2, 2, length = 100)
f_FS_α1 = Float64[]
f_FS_α01 = Float64[]
f_FS_α001 = Float64[]

for Kn in Kn_range
    prob = NonlinearProblem(compiled, Dict(compiled.Kn => Kn, compiled.α => 1.0))
    push!(f_FS_α1, solve(prob)[compiled.f_FS])

    prob = NonlinearProblem(compiled, Dict(compiled.Kn => Kn, compiled.α => 0.1))
    push!(f_FS_α01, solve(prob)[compiled.f_FS])

    prob = NonlinearProblem(compiled, Dict(compiled.Kn => Kn, compiled.α => 0.01))
    push!(f_FS_α001, solve(prob)[compiled.f_FS])
end

plot(Kn_range, f_FS_α1, label = "α = 1.0", xscale = :log10,
    xlabel = "Knudsen Number (Kn)", ylabel = "J/Jc (Correction Factor)",
    title = "Fuchs-Sutugin Transition Regime Correction")
plot!(Kn_range, f_FS_α01, label = "α = 0.1")
plot!(Kn_range, f_FS_α001, label = "α = 0.01")
savefig("fuchs_sutugin.svg");
nothing # hide
```

![Fuchs-Sutugin correction factor](fuchs_sutugin.svg)

The correction factor approaches 1 in the continuum limit (Kn → 0) and decreases in the kinetic regime (Kn → ∞). Lower accommodation coefficients (α) reduce the mass transfer rate.

### State Variables

```@example mass_transfer
using DataFrames, Symbolics, DynamicQuantities

sys = MassTransfer()
vars = ModelingToolkit.unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(DynamicQuantities.dimension(ModelingToolkit.get_unit(v)))
               for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example mass_transfer
params = ModelingToolkit.parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(DynamicQuantities.dimension(ModelingToolkit.get_unit(p)))
               for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example mass_transfer
equations(sys)
```

## Characteristic Timescales

The module also provides components for calculating characteristic timescales of various mass transport processes:

```@docs
GasDiffusionTimescale
AqueousDiffusionTimescale
InterfacialTimescale
ReactionTimescale
SolidEquilibrationTimescale
AqueousEquilibrationTimescale
```

### Timescale Comparison

```@example mass_transfer
using Aerosol

# Calculate timescales for a 10 μm droplet
sys_dg = GasDiffusionTimescale()
sys_da = AqueousDiffusionTimescale()

compiled_dg = mtkcompile(sys_dg)
compiled_da = mtkcompile(sys_da)

R_p = 1.0e-5  # 10 μm
D_g = 2.0e-5  # m²/s (gas-phase diffusivity)
D_aq = 1.0e-9  # m²/s (aqueous-phase diffusivity)

prob_dg = NonlinearProblem(compiled_dg, Dict(compiled_dg.R_p => R_p, compiled_dg.D_g => D_g))
τ_dg = solve(prob_dg)[compiled_dg.τ_dg]

prob_da = NonlinearProblem(compiled_da, Dict(compiled_da.R_p => R_p, compiled_da.D_aq => D_aq))
τ_da = solve(prob_da)[compiled_da.τ_da]

println("Gas-phase diffusion timescale: $(round(τ_dg * 1e6, digits=2)) μs")
println("Aqueous-phase diffusion timescale: $(round(τ_da * 1e3, digits=2)) ms")
println("Ratio τ_da/τ_dg: $(round(τ_da/τ_dg, digits=0))")
```

## Aqueous-Phase Chemistry Coupling

For droplet chemistry applications, the module provides components for mass transport limitations and coupled gas-aqueous dynamics:

```@docs
AqueousDiffusionReaction
MassTransportLimitation
DropletMassBalance
```

### Mass Transport Limitation Criteria

The following example calculates the limiting reaction rates below which mass transport does not limit aqueous-phase chemistry:

```@example mass_transfer
using Aerosol

sys = MassTransportLimitation()
compiled = mtkcompile(sys)

# Parameters for a 10 μm droplet at 298 K
prob = NonlinearProblem(compiled,
    Dict(
        compiled.R_p => 1.0e-5, compiled.D_g => 2.0e-5, compiled.D_aq => 1.0e-9,
        compiled.T => 298.15, compiled.α => 1.0, compiled.M_A => 0.029))
sol = solve(prob)

println("Gas-phase diffusion limit for k₁H*: $(sol[compiled.k1H_gas_limit]) mol/(m³·Pa·s)")
println("Aqueous-phase diffusion limit for k₁: $(sol[compiled.k1_aq_limit]) s⁻¹")
println("Interfacial transport limit for k₁H*: $(sol[compiled.k1H_interface_limit]) mol/(m³·Pa·s)")
```

## Analysis

### Comparison of Transition Regime Formulas

Figure 12.2 in Seinfeld & Pandis compares different transition regime formulas. Here we compare Fuchs-Sutugin and Dahneke:

```@example mass_transfer
sys_fs = FuchsSutugin()
sys_d = Dahneke()
compiled_fs = mtkcompile(sys_fs)
compiled_d = mtkcompile(sys_d)

Kn_range = 10 .^ range(-2, 2, length = 100)
f_FS = Float64[]
f_D = Float64[]

for Kn in Kn_range
    prob_fs = NonlinearProblem(compiled_fs, Dict(compiled_fs.Kn => Kn, compiled_fs.α => 1.0))
    push!(f_FS, solve(prob_fs)[compiled_fs.f_FS])

    prob_d = NonlinearProblem(compiled_d, Dict(compiled_d.Kn => Kn, compiled_d.α => 1.0))
    push!(f_D, solve(prob_d)[compiled_d.f_D])
end

plot(Kn_range, f_FS, label = "Fuchs-Sutugin", xscale = :log10,
    xlabel = "Knudsen Number (Kn)", ylabel = "J/Jc",
    title = "Comparison of Transition Regime Formulas (α = 1)")
plot!(Kn_range, f_D, label = "Dahneke", linestyle = :dash)
savefig("transition_comparison.svg");
nothing # hide
```

![Transition regime formula comparison](transition_comparison.svg)

Both formulas agree in the limiting regimes but show slight differences in the transition regime (Kn ~ 1).

### Particle Size Dependence

The mass transfer regime depends strongly on particle size through the Knudsen number:

```@example mass_transfer
sys = MassTransfer()
compiled = mtkcompile(sys)

# Calculate for particle diameters from 10 nm to 10 μm
D_p_range = 10 .^ range(-8, -5, length = 50)  # meters
R_p_range = D_p_range ./ 2

Kn_vals = Float64[]
f_FS_vals = Float64[]

for R_p in R_p_range
    prob = NonlinearProblem(compiled,
        Dict(
            compiled.R_p => R_p, compiled.α => 1.0,
            compiled.c_inf => 1.0e-6, compiled.c_s => 0.0))
    sol = solve(prob)
    push!(Kn_vals, sol[compiled.Kn])
    push!(f_FS_vals, sol[compiled.f_FS])
end

plot(D_p_range .* 1e6, Kn_vals, label = "Knudsen Number", xscale = :log10, yscale = :log10,
    xlabel = "Particle Diameter (μm)", ylabel = "Kn or J/Jc",
    title = "Particle Size Dependence")
plot!(D_p_range .* 1e6, f_FS_vals, label = "J/Jc (Fuchs-Sutugin)")
hline!([1.0], label = "Kn = 1 (transition)", linestyle = :dash, color = :gray)
savefig("size_dependence.svg");
nothing # hide
```

![Particle size dependence](size_dependence.svg)

Small particles (< 100 nm) are in the kinetic or transition regime where mass transfer is significantly reduced compared to the continuum limit.
