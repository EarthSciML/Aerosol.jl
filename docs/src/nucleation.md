# Nucleation

## Overview

Nucleation is the process by which gas-phase molecules spontaneously form small clusters (nuclei) that can grow into aerosol particles. This module implements classical nucleation theory for single-component (homogeneous) nucleation following the treatment in Seinfeld & Pandis (2006), Chapter 11.

The key outputs of the nucleation model are:

  - **Critical cluster size** (i*): The minimum number of molecules required for a stable cluster
  - **Critical radius** (r*): The radius of the critical cluster
  - **Nucleation rate** (J): The rate at which new particles form per unit volume per unit time
  - **Free energy barrier** (ΔG*): The energy barrier that must be overcome for nucleation to occur

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 2nd Edition, Chapter 11, pp. 489-536.

```@docs
Nucleation
CriticalCluster
ClassicalNucleationRate
WaterProperties
```

## Implementation

The implementation follows classical nucleation theory, which treats the formation of molecular clusters using thermodynamic arguments. The key assumption is the **capillarity approximation**: clusters of all sizes are assumed to have the same surface tension as bulk liquid.

### State Variables

```@example nucleation
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Aerosol
using OrdinaryDiffEq
using Plots

sys = Nucleation()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [ModelingToolkit.get_unit(v) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example nucleation
pars = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in pars],
    :Units => [ModelingToolkit.get_unit(p) for p in pars],
    :Description => [ModelingToolkit.getdescription(p) for p in pars]
)
```

### Equations

The system implements the following equations from Seinfeld & Pandis (2006):

```@example nucleation
eqs = equations(sys)
```

### Core Equations

The **saturation ratio** (Eq. 11.1) is the ratio of partial pressure to saturation vapor pressure:

```math
S = \frac{p_A}{p_A^s(T)}
```

The **critical cluster size** (Eq. 11.35) is determined by the competition between bulk and surface terms:

```math
i^* = \left(\frac{2\theta}{3\ln S}\right)^3
```

where θ is the dimensionless surface tension parameter (Eq. 11.27):

```math
\theta = \frac{(36\pi)^{1/3} v_1^{2/3} \sigma}{k_B T}
```

The **critical radius** (Eq. 11.52) follows from the Kelvin equation:

```math
r^* = \frac{2\sigma v_1}{k_B T \ln S}
```

The **free energy barrier** (Eq. 11.53) is:

```math
\Delta G^* = \frac{16\pi}{3} \frac{v_1^2 \sigma^3}{(k_B T \ln S)^2}
```

The **classical nucleation rate** (Eq. 11.47) gives the rate of new particle formation:

```math
J = \sqrt{\frac{2\sigma}{\pi m_1}} \frac{v_1 N_1^2}{S} \exp\left(-\frac{16\pi}{3} \frac{v_1^2 \sigma^3}{(k_B T)^3 (\ln S)^2}\right)
```

## Analysis

### Validation Against Table 11.1: Critical Cluster Properties for Water

The following analysis reproduces Table 11.1 from Seinfeld & Pandis (2006), showing critical cluster properties for water at 273 K and 298 K.

```@example nucleation
using OrdinaryDiffEq

# Water properties from Table 11.1
σ_273 = 75.6e-3  # N/m at 273 K
σ_298 = 72.0e-3  # N/m at 298 K
v_1 = 2.99e-29   # m³/molecule

sys = CriticalCluster()
compiled = mtkcompile(sys)

# Calculate for 273 K
results_273 = []
for S in [2, 3, 4, 5]
    prob = ODEProblem(compiled,
        Dict(
            compiled.T => 273.0,
            compiled.S => Float64(S),
            compiled.σ => σ_273,
            compiled.v_1 => v_1
        ),
        (0.0, 1.0))
    sol = solve(prob, Tsit5())
    r_star_angstrom = sol[compiled.r_star][1] * 1e10  # Convert to Å
    i_star = sol[compiled.i_star][1]
    push!(results_273, (
        S = S, r_star = round(r_star_angstrom, digits = 1), i_star = round(Int, i_star)))
end

# Calculate for 298 K
results_298 = []
for S in [2, 3, 4, 5]
    prob = ODEProblem(compiled,
        Dict(
            compiled.T => 298.0,
            compiled.S => Float64(S),
            compiled.σ => σ_298,
            compiled.v_1 => v_1
        ),
        (0.0, 1.0))
    sol = solve(prob, Tsit5())
    r_star_angstrom = sol[compiled.r_star][1] * 1e10
    i_star = sol[compiled.i_star][1]
    push!(results_298, (
        S = S, r_star = round(r_star_angstrom, digits = 1), i_star = round(Int, i_star)))
end

# Display results as a table
println("Table 11.1: Critical Number and Radius for Water Droplets")
println("="^70)
println("        T = 273 K                    T = 298 K")
println("S    r* (Å)    i*                r* (Å)    i*")
println("-"^70)
for (r273, r298) in zip(results_273, results_298)
    println("$(r273.S)    $(r273.r_star)       $(r273.i_star)                 $(r298.r_star)       $(r298.i_star)")
end

# Expected values from Table 11.1
println("\nExpected values from Seinfeld & Pandis Table 11.1:")
println("T=273K: S=2: r*=17.3Å, i*=726; S=3: r*=10.9Å, i*=182; S=4: r*=8.7Å, i*=87; S=5: r*=7.5Å, i*=58")
println("T=298K: S=2: r*=15.1Å, i*=482; S=3: r*=9.5Å, i*=121; S=4: r*=7.6Å, i*=60; S=5: r*=6.5Å, i*=39")
```

### Validation Against Table 11.4: Nucleation Rates for Water

This reproduces Table 11.4 from Seinfeld & Pandis (2006), showing nucleation rates for water at 293 K.

```@example nucleation
# Water properties at 293 K from Table 11.4
σ_293 = 72.75e-3      # N/m
v_1 = 2.99e-29        # m³/molecule
m_1 = 2.99e-26        # kg/molecule
p_sat_293 = 2336.5    # Pa

sys = Nucleation()
compiled = mtkcompile(sys)

println("Table 11.4: Homogeneous Nucleation Rate and Critical Cluster Size for Water at T = 293 K")
println("="^80)
println("S      i*      J (cm⁻³s⁻¹)           Expected J")
println("-"^80)

expected_J = Dict(2 => 5.02e-54, 3 => 1.76e-6, 4 => 1.05e6, 5 => 1.57e11,
    6 => 1.24e14, 7 => 8.99e15, 8 => 1.79e17, 9 => 1.65e18, 10 => 9.17e18)

for S in [2, 3, 4, 5, 6, 7, 8, 9, 10]
    prob = ODEProblem(compiled,
        Dict(
            compiled.T => 293.0,
            compiled.p_A => Float64(S) * p_sat_293,
            compiled.p_A_s => p_sat_293,
            compiled.σ => σ_293,
            compiled.v_1 => v_1,
            compiled.m_1 => m_1
        ),
        (0.0, 1.0))
    sol = solve(prob, Tsit5())
    J = sol[compiled.J][1] * 1e-6  # Convert from m⁻³s⁻¹ to cm⁻³s⁻¹
    i_star = sol[compiled.i_star][1]
    println("$(S)      $(round(Int, i_star))      $(round(J, sigdigits=3))           $(expected_J[S])")
end
```

### Figure: Critical Cluster Size vs. Saturation Ratio

This figure shows how the critical cluster size decreases with increasing saturation ratio, demonstrating why nucleation "switches on" sharply at a threshold supersaturation.

```@example nucleation
using Plots

S_range = 1.5:0.1:10.0
i_stars = Float64[]
r_stars = Float64[]

sys = CriticalCluster()
compiled = mtkcompile(sys)

for S in S_range
    prob = ODEProblem(compiled,
        Dict(
            compiled.T => 293.0,
            compiled.S => S,
            compiled.σ => 72.75e-3,
            compiled.v_1 => 2.99e-29
        ),
        (0.0, 1.0))
    sol = solve(prob, Tsit5())
    push!(i_stars, sol[compiled.i_star][1])
    push!(r_stars, sol[compiled.r_star][1] * 1e10)  # Convert to Å
end

p1 = plot(
    S_range, i_stars, xlabel = "Saturation Ratio S", ylabel = "Critical Cluster Size i*",
    label = "i*", lw = 2, legend = :topright, yscale = :log10,
    title = "Critical Cluster Size vs. Saturation Ratio (T=293K)")

p2 = plot(
    S_range, r_stars, xlabel = "Saturation Ratio S", ylabel = "Critical Radius r* (Å)",
    label = "r*", lw = 2, legend = :topright,
    title = "Critical Radius vs. Saturation Ratio (T=293K)")

plot(p1, p2, layout = (1, 2), size = (900, 400))
```

### Figure: Nucleation Rate vs. Saturation Ratio (analogous to Table 11.4)

This figure illustrates the extreme sensitivity of nucleation rate to saturation ratio.

```@example nucleation
S_range = 2.0:0.5:10.0
J_rates = Float64[]

sys = Nucleation()
compiled = mtkcompile(sys)

for S in S_range
    prob = ODEProblem(compiled,
        Dict(
            compiled.T => 293.0,
            compiled.p_A => S * 2336.5,
            compiled.p_A_s => 2336.5,
            compiled.σ => 72.75e-3,
            compiled.v_1 => 2.99e-29,
            compiled.m_1 => 2.99e-26
        ),
        (0.0, 1.0))
    sol = solve(prob, Tsit5())
    J = sol[compiled.J][1] * 1e-6  # Convert to cm⁻³s⁻¹
    push!(J_rates, max(J, 1e-60))  # Floor at 1e-60 for log scale
end

plot(
    S_range, J_rates, xlabel = "Saturation Ratio S", ylabel = "Nucleation Rate J (cm⁻³s⁻¹)",
    label = "J", lw = 2, yscale = :log10,
    title = "Classical Nucleation Rate vs. Saturation Ratio (Water at 293K)",
    ylims = (1e-60, 1e25), legend = :bottomright)
```

### Temperature Dependence of Critical Properties

This analysis shows how temperature affects critical cluster properties at a fixed saturation ratio.

```@example nucleation
T_range = 250.0:5.0:320.0
S = 4.0

i_stars_T = Float64[]
r_stars_T = Float64[]

# Approximate temperature-dependent surface tension: σ ≈ 75.6 - 0.1454*(T-273) mN/m
σ_func(T) = (75.6 - 0.1454 * (T - 273)) * 1e-3

sys = CriticalCluster()
compiled = mtkcompile(sys)

for T in T_range
    prob = ODEProblem(compiled,
        Dict(
            compiled.T => T,
            compiled.S => S,
            compiled.σ => σ_func(T),
            compiled.v_1 => 2.99e-29
        ),
        (0.0, 1.0))
    sol = solve(prob, Tsit5())
    push!(i_stars_T, sol[compiled.i_star][1])
    push!(r_stars_T, sol[compiled.r_star][1] * 1e10)
end

p1 = plot(
    T_range, i_stars_T, xlabel = "Temperature (K)", ylabel = "Critical Cluster Size i*",
    label = "i*", lw = 2, title = "Critical Cluster Size vs. Temperature (S=4)")

p2 = plot(
    T_range, r_stars_T, xlabel = "Temperature (K)", ylabel = "Critical Radius r* (Å)",
    label = "r*", lw = 2, title = "Critical Radius vs. Temperature (S=4)")

plot(p1, p2, layout = (1, 2), size = (900, 400))
```

## Key Behaviors

### Saturation Ratio Dependence

The nucleation rate is extremely sensitive to the saturation ratio S:

  - At S < 2, nucleation is essentially zero for atmospheric conditions
  - The transition from negligible to significant nucleation occurs over a narrow range of S
  - At high S (> 5), the rate increases by many orders of magnitude per unit increase in S

### Temperature Dependence

Temperature affects nucleation through multiple pathways:

  - Higher T reduces the free energy barrier, increasing J
  - Higher T increases vapor pressure, affecting S
  - Surface tension typically decreases with T

### Critical Cluster Size

The critical cluster size i* decreases with increasing S:

  - At low S, very large clusters are needed before they become stable
  - At high S, even small clusters can grow spontaneously
  - This explains why nucleation "switches on" sharply at a threshold S

## Limitations

 1. **Capillarity approximation**: The model assumes clusters of all sizes have bulk liquid surface tension, which becomes inaccurate for small clusters (i* < 20-30).

 2. **Single-component only**: This implementation covers homogeneous single-component nucleation. Binary nucleation (e.g., H₂SO₄-H₂O) requires additional equations not yet implemented.
 3. **Steady-state assumption**: The classical nucleation rate assumes the cluster distribution has reached steady state.
 4. **No kinetic effects**: The model does not account for non-equilibrium effects or condensation/evaporation kinetics beyond the equilibrium rate constants.
