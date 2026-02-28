# Stochastic Collection Equation

## Overview

This module implements the two-moment method for solving the stochastic collection equation (SCE) for cloud droplet coalescence. The SCE describes the evolution of a drop size distribution due to gravitational collection, where larger drops fall faster and sweep up smaller ones.

The method conserves two moments per mass category: number concentration ``N_k`` and mass concentration ``M_k``. Closure is achieved using a nondimensional parameter ``\bar{\xi}_p`` that relates higher-order moments to these two tracked moments, avoiding the need for weighting functions. The approach is more accurate and computationally efficient than single-moment methods such as Bleck's (1970) algorithm.

**Reference**: Tzivion, S., Feingold, G., and Levin, Z. (1989) "The Evolution of Raindrop Spectra. Part II: Collisional Collection/Breakup and Evaporation in a Rainshaft", *Journal of the Atmospheric Sciences*, 46(21), 3312-3327.

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

### Constant Kernel — Fractional Mass Distribution (cf. Figure 1a)

For the constant kernel ``K(x,y) = K_0 = 10^{-4}`` cm³ s⁻¹ (= ``10^{-10}`` m³ s⁻¹), with an exponential initial distribution (``N_0 = 300`` cm⁻³, LWC = 1 g m⁻³), the following figure shows the fractional mass distribution ``M_k / \text{LWC}`` at ``t = 0`` and ``t = 60`` minutes of collection. The spectrum shifts to larger categories as coalescence proceeds, consistent with Figure 1a in the paper.

```@example sce
using OrdinaryDiffEq
using Plots

I = 8; x1 = 1.6e-14; K0_val = 1e-10
sys_const = StochasticCollectionCoalescence(; I=I, kernel_type=:constant, kernel_params=Dict(:K0 => K0_val))
compiled_const = mtkcompile(sys_const)

N0 = 300e6; LWC = 1e-3
xbar = LWC / N0
u0 = Dict{Any,Float64}()
for k in 1:I
    xk = x1 * 2.0^(k - 1); xk1 = x1 * 2.0^k
    u0[compiled_const.Nk[k]] = max(N0 * (exp(-xk / xbar) - exp(-xk1 / xbar)), 1e-20)
    u0[compiled_const.Mk[k]] = max(N0 * xbar * ((xk / xbar + 1) * exp(-xk / xbar) - (xk1 / xbar + 1) * exp(-xk1 / xbar)), 1e-30)
end

prob = ODEProblem(compiled_const, u0, (0.0, 3600.0))
sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-12)

# Plot fractional mass at t=0 and t=60 min
frac_mass_0 = [sol[compiled_const.Mk[k]][1] / LWC for k in 1:I]
frac_mass_60 = [sol[compiled_const.Mk[k]][end] / LWC for k in 1:I]

p = plot(1:I, frac_mass_0, label="0 min", xlabel="Category (k)",
    ylabel="Fractional Mass (Mₖ/LWC)", title="Constant Kernel: Fractional Mass Distribution",
    marker=:circle, lw=2)
plot!(p, 1:I, frac_mass_60, label="60 min", marker=:square, lw=2)
p
```

### Constant Kernel — Number Concentration (cf. Figure 1b)

Category number concentrations ``N_k`` at ``t = 0`` and ``t = 60`` min of collection. As coalescence proceeds, small-drop categories are depleted and larger categories gain particles.

```@example sce
Nk_0 = [sol[compiled_const.Nk[k]][1] for k in 1:I]
Nk_60 = [sol[compiled_const.Nk[k]][end] for k in 1:I]

p = plot(1:I, Nk_0, label="0 min", xlabel="Category (k)",
    ylabel="Nₖ (m⁻³)", title="Constant Kernel: Number Concentration",
    marker=:circle, lw=2, yscale=:log10)
plot!(p, 1:I, Nk_60, label="60 min", marker=:square, lw=2)
p
```

### Golovin Kernel — Fractional Mass Distribution (cf. Figure 2a)

For the Golovin kernel ``K(x,y) = C(x+y)`` with ``C = 1500`` cm³ s⁻¹ g⁻¹ (= ``1.5 \times 10^{-3}`` m³ s⁻¹ kg⁻¹), solutions are shown at ``t = 0`` and ``t = 30`` minutes.

```@example sce
C_val = 1.5e-3
sys_gol = StochasticCollectionCoalescence(; I=I, kernel_type=:golovin, kernel_params=Dict(:C => C_val))
compiled_gol = mtkcompile(sys_gol)

u0_gol = Dict{Any,Float64}()
for k in 1:I
    xk = x1 * 2.0^(k - 1); xk1 = x1 * 2.0^k
    u0_gol[compiled_gol.Nk[k]] = max(N0 * (exp(-xk / xbar) - exp(-xk1 / xbar)), 1e-20)
    u0_gol[compiled_gol.Mk[k]] = max(N0 * xbar * ((xk / xbar + 1) * exp(-xk / xbar) - (xk1 / xbar + 1) * exp(-xk1 / xbar)), 1e-30)
end

prob_gol = ODEProblem(compiled_gol, u0_gol, (0.0, 1800.0))
sol_gol = solve(prob_gol, Tsit5(); reltol=1e-8, abstol=1e-12)

frac_mass_gol_0 = [sol_gol[compiled_gol.Mk[k]][1] / LWC for k in 1:I]
frac_mass_gol_30 = [sol_gol[compiled_gol.Mk[k]][end] / LWC for k in 1:I]

p = plot(1:I, frac_mass_gol_0, label="0 min", xlabel="Category (k)",
    ylabel="Fractional Mass (Mₖ/LWC)", title="Golovin Kernel: Fractional Mass Distribution",
    marker=:circle, lw=2)
plot!(p, 1:I, frac_mass_gol_30, label="30 min", marker=:square, lw=2)
p
```

### Golovin Kernel — Number Concentration (cf. Figure 2b)

Category number concentrations ``N_k`` for the Golovin kernel at ``t = 0`` and ``t = 30`` min of collection, showing how coalescence redistributes particles among size categories.

```@example sce
Nk_gol_0 = [sol_gol[compiled_gol.Nk[k]][1] for k in 1:I]
Nk_gol_30 = [sol_gol[compiled_gol.Nk[k]][end] for k in 1:I]

p = plot(1:I, Nk_gol_0, label="0 min", xlabel="Category (k)",
    ylabel="Nₖ (m⁻³)", title="Golovin Kernel: Number Concentration",
    marker=:circle, lw=2, yscale=:log10)
plot!(p, 1:I, Nk_gol_30, label="30 min", marker=:square, lw=2)
p
```

### Mass Conservation

A key advantage of the two-moment method is conservation of both number and mass. The following shows mass conservation over time for both kernels:

```@example sce
times_const = sol.t
Mk_const = [sol[compiled_const.Mk[k]] for k in 1:I]
M_total_const = [sum(Mk_const[k][j] for k in 1:I) for j in eachindex(times_const)]

times_gol = sol_gol.t
Mk_gol = [sol_gol[compiled_gol.Mk[k]] for k in 1:I]
M_total_gol = [sum(Mk_gol[k][j] for k in 1:I) for j in eachindex(times_gol)]

p = plot(times_const ./ 60, M_total_const ./ LWC,
    label="Constant kernel", xlabel="Time (min)",
    ylabel="M_total / LWC₀", title="Mass Conservation",
    lw=2, ylims=(0.0, 1.2))
plot!(p, times_gol ./ 60, M_total_gol ./ LWC,
    label="Golovin kernel", lw=2, ls=:dash)
p
```
