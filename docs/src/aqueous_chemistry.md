# Aqueous Chemistry

## Overview

Atmospheric aqueous-phase chemistry module implementing equations from Seinfeld & Pandis Chapter 7: "Atmospheric Aqueous-Phase Chemistry".

Aqueous-phase chemistry occurs in cloud droplets, fog droplets, and deliquesced aerosol particles. These processes are critical for:

  - Formation of secondary aerosol components (especially sulfate)
  - pH regulation in atmospheric water
  - Redistribution of species between gas and aqueous phases
  - Oxidation of dissolved gases

This module provides ModelingToolkit.jl systems for:

  - Henry's law gas-liquid equilibrium
  - Aqueous dissociation equilibria (CO2, SO2, NH3, HNO3, H2O2, O3)
  - S(IV) to S(VI) oxidation kinetics
  - Combined cloud droplet chemistry

**Reference**: Seinfeld, J.H. and Pandis, S.N., "Atmospheric Chemistry and Physics: From Air Pollution to Climate Change", 3rd Edition, Wiley, 2016, Chapter 7.

```@docs
CloudChemistry
```

## Implementation

All pressures use SI units (Pa). Henry's law constants from the textbook (M/atm) are internally converted to M/Pa by dividing by 101325. Concentrations are in mol/L, temperatures in K.

### Henry's Law Components

Henry's law describes the equilibrium partitioning of gases between the gas phase and aqueous phase.

```@docs
HenrysLaw
HenrysLawTemperature
EffectiveHenrysLaw
```

### Aqueous Equilibria

```@docs
WaterEquilibrium
CO2Equilibria
SO2Equilibria
NH3Equilibria
HNO3Equilibria
H2O2Equilibria
O3Equilibria
AqueousEquilibria
```

### Sulfate Formation (S(IV) Oxidation)

```@docs
SulfateFormationO3
SulfateFormationH2O2
SulfateFormationFe
SulfateFormationMn
SulfateFormationFeMn
SulfateFormation
```

### Cloud Chemistry Systems

```@docs
CloudChemistryFixedpH
CloudChemistryODE
```

### Utility Functions

```@docs
effective_henrys_constant
aqueous_fraction
distribution_factor
henrys_constant_at_T
rate_to_ppb_hr
rate_to_percent_hr
so2_lifetime
```

### State Variables

```@example aqueous_chem
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Aerosol

sys = CloudChemistry()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example aqueous_chem
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example aqueous_chem
eqs = equations(sys)
```

## Analysis

### S(IV) Mole Fraction Diagram

The distribution of S(IV) species (SO₂·H₂O, HSO₃⁻, SO₃²⁻) depends on pH.
The mole fractions (alpha values) are computed from the equilibrium constants
(S&P Eq. 7.43-7.45).

```@example aqueous_chem
using Plots

K_s1 = Aerosol.K_S1_298
K_s2 = Aerosol.K_S2_298

pH_range = range(0, 10, length = 200)
alpha_0 = zeros(length(pH_range))
alpha_1 = zeros(length(pH_range))
alpha_2 = zeros(length(pH_range))

for (i, pH) in enumerate(pH_range)
    H_plus = 10.0^(-pH)
    denom = H_plus^2 + K_s1 * H_plus + K_s1 * K_s2
    alpha_0[i] = H_plus^2 / denom
    alpha_1[i] = K_s1 * H_plus / denom
    alpha_2[i] = K_s1 * K_s2 / denom
end

plot(pH_range, alpha_0, label = "α₀ (SO₂·H₂O)", lw = 2)
plot!(pH_range, alpha_1, label = "α₁ (HSO₃⁻)", lw = 2)
plot!(pH_range, alpha_2, label = "α₂ (SO₃²⁻)", lw = 2)
xlabel!("pH")
ylabel!("Mole fraction")
title!("S(IV) Species Distribution (S&P Fig. 7.7)")
```

### Effective Henry's Law Constant vs pH

The effective Henry's law constant H* for SO₂ increases with pH due to
dissociation of dissolved SO₂ (S&P Eq. 7.41).

```@example aqueous_chem
ATM = Aerosol.ATM_TO_PA
H_SO2_298 = 1.23 / ATM  # M/Pa

pH_range = range(0, 8, length = 200)
H_eff_vals = zeros(length(pH_range))
H_eff_atm = zeros(length(pH_range))

for (i, pH) in enumerate(pH_range)
    H_plus = 10.0^(-pH)
    H_eff_vals[i] = effective_henrys_constant(H_SO2_298, K_s1, K_s2, H_plus)
    H_eff_atm[i] = H_eff_vals[i] * ATM  # Convert back to M/atm for display
end

plot(pH_range, H_eff_atm, lw = 2, yscale = :log10, label = "H* SO₂")
hline!([1.23], ls = :dash, label = "H (intrinsic)", color = :gray)
xlabel!("pH")
ylabel!("H* (M/atm)")
title!("Effective Henry's Law Constant for SO₂ (S&P Eq. 7.41)")
```

### S(IV) Oxidation Rates vs pH

Comparison of the three main S(IV) oxidation pathways as a function of pH.
This figure illustrates the pH-dependent dominance of different oxidants
(cf. S&P Section 7.5).

```@example aqueous_chem
ATM = Aerosol.ATM_TO_PA
H_SO2 = 1.23 / ATM
H_O3 = 1.1e-2 / ATM
H_H2O2 = 1.0e5 / ATM
K_s1 = Aerosol.K_S1_298
K_s2 = Aerosol.K_S2_298

# Typical atmospheric conditions
p_SO2 = 1e-9 * ATM  # 1 ppb
p_O3 = 50e-9 * ATM  # 50 ppb
p_H2O2 = 1e-9 * ATM # 1 ppb

pH_range = range(1, 7, length = 200)
R_o3 = zeros(length(pH_range))
R_h2o2 = zeros(length(pH_range))

k0 = Aerosol.K0_O3
k1 = Aerosol.K1_O3
k2 = Aerosol.K2_O3
k_h2o2 = 7.5e7
K_eq = 13.0

for (i, pH) in enumerate(pH_range)
    H_plus = 10.0^(-pH)
    SO2_aq = H_SO2 * p_SO2
    HSO3_m = K_s1 * SO2_aq / H_plus
    SO3_2m = K_s2 * HSO3_m / H_plus
    O3_aq = H_O3 * p_O3
    H2O2_aq = H_H2O2 * p_H2O2

    R_o3[i] = (k0 * SO2_aq + k1 * HSO3_m + k2 * SO3_2m) * O3_aq
    R_h2o2[i] = k_h2o2 * H_plus * H2O2_aq * HSO3_m / (1 + K_eq * H_plus)
end

plot(pH_range, R_o3, label = "O₃ pathway", lw = 2, yscale = :log10)
plot!(pH_range, R_h2o2, label = "H₂O₂ pathway", lw = 2)
xlabel!("pH")
ylabel!("Rate (mol/L/s)")
title!("S(IV) Oxidation Rates vs pH (S&P Section 7.5)")
```

### Temperature Dependence of Henry's Law Constants

Henry's law constants increase at lower temperatures for all species
(exothermic dissolution), following the van't Hoff equation (S&P Eq. 7.5).

```@example aqueous_chem
ATM = Aerosol.ATM_TO_PA
T_range = range(273, 313, length = 100)
species_list = [:SO2, :CO2, :NH3, :H2O2, :O3]
labels = ["SO₂", "CO₂", "NH₃", "H₂O₂", "O₃"]

plt = plot(yscale = :log10)
for (sp, lab) in zip(species_list, labels)
    H_298 = Aerosol.HENRY_CONSTANTS_298[sp]
    dH = Aerosol.DELTA_H_DISSOLUTION[sp]
    H_vals = [henrys_constant_at_T(H_298, dH, T) * ATM for T in T_range]
    plot!(plt, collect(T_range), H_vals, label = lab, lw = 2)
end
xlabel!(plt, "Temperature (K)")
ylabel!(plt, "H (M/atm)")
title!(plt, "Henry's Law Constants vs Temperature (S&P Eq. 7.5)")
plt
```
