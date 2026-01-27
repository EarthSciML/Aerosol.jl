# Aqueous Chemistry

Atmospheric aqueous-phase chemistry module implementing equations from Seinfeld & Pandis Chapter 7: "Atmospheric Aqueous-Phase Chemistry".

This module provides ModelingToolkit.jl systems for:
- Henry's law gas-liquid equilibrium
- Aqueous dissociation equilibria (CO2, SO2, NH3, HNO3, H2O2, O3)
- S(IV) to S(VI) oxidation kinetics
- Combined cloud droplet chemistry

## Overview

Aqueous-phase chemistry occurs in cloud droplets, fog droplets, and deliquesced aerosol particles. These processes are critical for:
- Formation of secondary aerosol components (especially sulfate)
- pH regulation in atmospheric water
- Redistribution of species between gas and aqueous phases
- Oxidation of dissolved gases

## Henry's Law Components

Henry's law describes the equilibrium partitioning of gases between the gas phase and aqueous phase.

### Basic Henry's Law

```@docs
HenrysLaw
```

Implements the fundamental relationship:
```math
[A(aq)] = H_A \cdot p_A
```

where `[A(aq)]` is the aqueous concentration (mol/L), `H_A` is Henry's law constant (M/atm), and `p_A` is the partial pressure (atm).

### Temperature-Dependent Henry's Law

```@docs
HenrysLawTemperature
```

Includes temperature dependence via the van't Hoff equation:
```math
H_A(T) = H_A(T_{ref}) \exp\left(\frac{\Delta H_{diss}}{R}\left(\frac{1}{T_{ref}} - \frac{1}{T}\right)\right)
```

Also calculates the distribution factor and aqueous fraction.

### Effective Henry's Law

```@docs
EffectiveHenrysLaw
```

For species that dissociate in water, the effective Henry's law constant accounts for all dissolved forms:
```math
H^* = H \left(1 + \frac{K_1}{[H^+]} + \frac{K_1 K_2}{[H^+]^2}\right)
```

## Aqueous Equilibria

### Water Autoionization

```@docs
WaterEquilibrium
```

The fundamental equilibrium:
```math
H_2O \rightleftharpoons H^+ + OH^-
```

### Carbon Dioxide System

```@docs
CO2Equilibria
```

The carbonate system:
```math
\begin{align}
CO_2 \cdot H_2O &\rightleftharpoons H^+ + HCO_3^- \\
HCO_3^- &\rightleftharpoons H^+ + CO_3^{2-}
\end{align}
```

### Sulfur Dioxide System

```@docs
SO2Equilibria
```

The S(IV) system:
```math
\begin{align}
SO_2 \cdot H_2O &\rightleftharpoons H^+ + HSO_3^- \\
HSO_3^- &\rightleftharpoons H^+ + SO_3^{2-}
\end{align}
```

### Other Species

```@docs
NH3Equilibria
HNO3Equilibria
H2O2Equilibria
O3Equilibria
```

### Combined Equilibria System

```@docs
AqueousEquilibria
```

## Sulfate Formation (S(IV) Oxidation)

Multiple pathways oxidize S(IV) to S(VI) (sulfate) in aqueous solution.

### Ozone Pathway

```@docs
SulfateFormationO3
```

Fast reaction, especially at high pH:
```math
R_{O_3} = (k_0[SO_2 \cdot H_2O] + k_1[HSO_3^-] + k_2[SO_3^{2-}])[O_3(aq)]
```

### Hydrogen Peroxide Pathway

```@docs
SulfateFormationH2O2
```

Dominant at low pH in clouds.

### Metal-Catalyzed Pathways

```@docs
SulfateFormationFe
SulfateFormationMn
SulfateFormationFeMn
```

Fe(III) and Mn(II) catalyze S(IV) oxidation by dissolved oxygen. The Fe/Mn combination shows synergistic enhancement.

### Combined Sulfate Formation

```@docs
SulfateFormation
```

Combines all oxidation pathways.

## Cloud Chemistry Systems

### Fixed pH Cloud Chemistry

```@docs
CloudChemistryFixedpH
```

Simplified system assuming constant pH.

### Full Cloud Chemistry System

```@docs
CloudChemistry
```

Complete system with pH calculation.

### ODE Formulation

```@docs
CloudChemistryODE
```

Time-dependent formulation for kinetic simulations.

## Utility Functions

```@docs
effective_henrys_constant
aqueous_fraction
distribution_factor
henrys_constant_at_T
rate_to_ppb_hr
rate_to_percent_hr
so2_lifetime
```

## Physical Constants

The module exports physical constants from Seinfeld & Pandis tables:

### Henry's Law Constants at 298 K
- `HENRY_CONSTANTS_298` - Dictionary of Henry's law constants (M/atm)
- `DELTA_H_DISSOLUTION` - Dictionary of heats of dissolution (J/mol)

### Equilibrium Constants at 298 K
- `K_W_298` - Water dissociation constant (1.0×10⁻¹⁴ M²)
- `K_C1_298`, `K_C2_298` - CO₂ system constants
- `K_S1_298`, `K_S2_298` - SO₂ system constants
- `K_A1_298` - NH₃ equilibrium constant

### Oxidation Rate Constants
- `K0_O3`, `K1_O3`, `K2_O3` - O₃ oxidation rate constants
- `K_H2O2` - H₂O₂ oxidation rate constant
- `K_FE_SIV`, `K_MN_SIV`, `K_FEMN_SIV` - Metal catalysis rate constants

## Units Convention

!!! note "Units Convention"
    Units are documented in descriptions but not enforced via DynamicQuantities due to non-SI units (atm) used in atmospheric chemistry conventions. All concentrations are in mol/L (M), pressures in atm, temperatures in K.

## Reference

Seinfeld, J.H. and Pandis, S.N., "Atmospheric Chemistry and Physics: From Air Pollution to Climate Change", 3rd Edition, Wiley, 2016, Chapter 7.
