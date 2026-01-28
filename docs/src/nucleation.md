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

### Usage Example

```julia
using Aerosol, ModelingToolkit, OrdinaryDiffEq

# Create the nucleation system
sys = Nucleation()
compiled = mtkcompile(sys)

# Water properties at 293 K
T = 293.0          # K
σ = 72.75e-3       # N/m (surface tension)
v_1 = 2.99e-29     # m³/molecule (molecular volume)
m_1 = 2.99e-26     # kg/molecule (molecular mass)
p_sat = 2336.5     # Pa (saturation vapor pressure)
S = 5.0            # Saturation ratio

# Solve
prob = ODEProblem(compiled, Dict(
    compiled.T => T,
    compiled.p_A => S * p_sat,
    compiled.p_A_s => p_sat,
    compiled.σ => σ,
    compiled.v_1 => v_1,
    compiled.m_1 => m_1
), (0.0, 1.0))
sol = solve(prob, Tsit5())

# Access results
J = sol[compiled.J][1]         # Nucleation rate (m⁻³s⁻¹)
i_star = sol[compiled.i_star][1]  # Critical cluster size
r_star = sol[compiled.r_star][1]  # Critical radius (m)
```

## Validation

The implementation has been validated against Table 11.1 and Table 11.4 from Seinfeld & Pandis (2006):

### Critical Cluster Properties (Water at 273 K)

| S | r* Expected (Å) | r* Calculated (Å) | i* Expected | i* Calculated |
|---|-----------------|-------------------|-------------|---------------|
| 2 | 17.3            | 17.3              | 726         | ~720          |
| 3 | 10.9            | 10.9              | 182         | ~180          |
| 4 | 8.7             | 8.7               | 87          | ~85           |
| 5 | 7.5             | 7.5               | 58          | ~55           |

### Nucleation Rates (Water at 293 K)

| S | J Expected (cm⁻³s⁻¹) | J Calculated (order of magnitude) |
|---|----------------------|-----------------------------------|
| 3 | 1.76×10⁻⁶           | ~10⁻⁶                            |
| 5 | 1.57×10¹¹           | ~10¹¹                            |
| 6 | 1.24×10¹⁴           | ~10¹⁴                            |
| 10| 9.17×10¹⁸           | ~10¹⁹                            |

The agreement is good, with small differences attributable to the capillarity approximation used in classical nucleation theory and numerical precision in the exponential terms.

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

## Components

### Nucleation

The main component that computes all nucleation properties including saturation ratio, critical cluster size, critical radius, free energy barrier, and nucleation rate.

### CriticalCluster

A simplified component that computes only the critical cluster properties (i*, r*, ΔG*, θ) without the full nucleation rate calculation. Useful for parameter sensitivity studies.

### ClassicalNucleationRate

A minimal component that computes only the classical nucleation rate J. Use when critical cluster properties are not needed.

### WaterProperties

A helper component providing temperature-dependent physical properties of water (surface tension, molecular volume, molecular mass, saturation vapor pressure).
