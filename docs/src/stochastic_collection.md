# Stochastic Collection Equation

## Overview

This module implements the two-moment method for solving the stochastic collection equation (SCE) for cloud droplet coalescence. The SCE describes the evolution of a drop size distribution due to gravitational collection, where larger drops fall faster and sweep up smaller ones.

The method conserves two moments per mass category: number concentration ``N_k`` and mass concentration ``M_k``. Closure is achieved using a nondimensional parameter ``\bar{\xi}_p`` that relates higher-order moments to these two tracked moments, avoiding the need for weighting functions. The approach is more accurate and computationally efficient than single-moment methods such as Bleck's (1970) algorithm.

**Reference**:
- Tzivion, S., Feingold, G., and Levin, Z. (1989) "The Evolution of Raindrop Spectra. Part II: Collisional Collection/Breakup and Evaporation in a Rainshaft", *Journal of the Atmospheric Sciences*, 46(21), 3312-3327.

## Implementation Scope and Limitations

This implementation covers only the **collisional collection** process described in the 1989 paper. The following processes from the paper are **not yet implemented**:

- **Breakup processes**: Critical for realistic evolution of larger raindrops (>2-3 mm diameter)
- **Evaporation processes**: Important for smaller droplets and subsaturated environments
- **Sedimentation processes**: Affects vertical distribution and spectral evolution with altitude

For a complete reproduction of the paper's results, these additional physics processes would need to be implemented and coupled with the collection equation.

```@docs
StochasticCollectionCoalescence
```

## Implementation

The continuous drop spectrum is divided into ``I`` discrete mass categories with geometrically increasing boundaries (``x_{k+1} = 2 x_k``, i.e. ``p = 2``). Two moments are tracked per category:
- ``N_k``: number concentration (m⁻³)
- ``M_k``: mass concentration (kg m⁻³)

The moment equations follow Eq. (9a,b) in the paper, with incomplete category integrals evaluated using a linear distribution function approximation (Eq. 11-13) and positivity constraints (Eq. 15a,b).

### State Variables

```@example sce
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Aerosol

sys = StochasticCollectionCoalescence(; I=5, kernel_type=:constant, kernel_params=Dict(:K0 => 1e-10))

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Equations

```@example sce
eqs = equations(sys)
```

## Analysis

### Constant Kernel — Spectral Evolution

For the constant kernel ``K(x,y) = K_0 = 10^{-10}`` m³ s⁻¹, with an exponential initial distribution (``N_0 = 300 \times 10^6`` m⁻³, LWC = 1 g m⁻³), the spectrum shifts to larger categories as coalescence proceeds, consistent with the paper's results.

```@example sce
using OrdinaryDiffEq
using Plots

I = 6; x1 = 1.6e-14; K0_val = 1e-10
sys_const = StochasticCollectionCoalescence(; I=I, kernel_type=:constant, kernel_params=Dict(:K0 => K0_val))
compiled_const = mtkcompile(sys_const)

N0 = 300e6; LWC = 1e-3; xbar = LWC / N0
u0 = Dict{Any,Float64}()
for k in 1:I
    xk = x1 * 2.0^(k - 1); xk1 = x1 * 2.0^k
    u0[compiled_const.Nk[k]] = max(N0 * (exp(-xk / xbar) - exp(-xk1 / xbar)), 1e-20)
    u0[compiled_const.Mk[k]] = max(N0 * xbar * ((xk / xbar + 1) * exp(-xk / xbar) - (xk1 / xbar + 1) * exp(-xk1 / xbar)), 1e-30)
end

prob = ODEProblem(compiled_const, u0, (0.0, 3600.0))
sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-12)

# Fractional mass distribution
frac_mass_0 = [sol[compiled_const.Mk[k]][1] / LWC for k in 1:I]
frac_mass_60 = [sol[compiled_const.Mk[k]][end] / LWC for k in 1:I]

plot(1:I, frac_mass_0, label="0 min", xlabel="Category (k)",
    ylabel="Fractional Mass", marker=:circle, lw=2, size=(400,250))
plot!(1:I, frac_mass_60, label="60 min", marker=:square, lw=2)
```

The total number concentration decreases from initial value due to coalescence, while total mass is approximately conserved (some mass may leak through the top category).

### Golovin Kernel

For the Golovin kernel ``K(x,y) = C(x+y)`` with ``C = 1.5 \times 10^{-3}`` m³ s⁻¹ kg⁻¹, the implementation supports mass-dependent collection rates that accelerate coalescence for larger droplets.

```@example sce
# Golovin kernel example (reduced output)
sys_gol = StochasticCollectionCoalescence(; I=4, kernel_type=:golovin, kernel_params=Dict(:C => 1.5e-3))
length(equations(sys_gol))  # Shows 8 equations (4 categories × 2 moments)
```

The two-moment method conserves both number concentration (decreases due to coalescence) and total mass (approximately conserved, with some loss through the top category boundary for finite category systems).
