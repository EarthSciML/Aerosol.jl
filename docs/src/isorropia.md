# ISORROPIA II

## Overview

ISORROPIA II is a computationally efficient thermodynamic equilibrium model for K⁺–Ca²⁺–Mg²⁺–NH₄⁺–Na⁺–SO₄²⁻–NO₃⁻–Cl⁻–H₂O aerosol systems. The model calculates the equilibrium partitioning between gas, liquid, and solid phases for multicomponent inorganic aerosols under given temperature, relative humidity, and total species concentrations.

**Reference**: Fountoukis, C. and Nenes, A., 2007. ISORROPIA II: a computationally efficient thermodynamic equilibrium model for K⁺–Ca²⁺–Mg²⁺–NH₄⁺–Na⁺–SO₄²⁻–NO₃⁻–Cl⁻–H₂O aerosols. *Atmospheric Chemistry and Physics*, 7(17), pp.4639-4659. [https://doi.org/10.5194/acp-7-4639-2007](https://doi.org/10.5194/acp-7-4639-2007)

```@docs
Isorropia
```

## Implementation

The ISORROPIA II model is implemented using ModelingToolkit.jl and includes:

- **Thermodynamic equilibrium constants** for 27 chemical reactions (Table 2 of reference paper)
- **Activity coefficient calculations** using multicomponent Debye-Hückel and Kusik-Meissner formulations
- **Aerosol type classification** based on sulfate ratio and crustal species content
- **Phase partitioning** between gas, aqueous, and solid phases
- **Temperature dependence** for all equilibrium constants

### Key Features

- Supports both forward and reverse problems
- Handles 5 aerosol types based on composition regimes
- Includes deliquescence relative humidity calculations
- Comprehensive test coverage against FORTRAN reference implementation

### Usage Example

```julia
using Aerosol
using ModelingToolkit
using OrdinaryDiffEq

# Create ISORROPIA system
@named isrpa = Isorropia()
sys = mtkcompile(isrpa)

# Set conditions (T=298K, RH=50%)
prob = ODEProblem(sys, [], (0.0, 1.0), [sys.T => 298.15, sys.RH => 0.5])
sol = solve(prob, Rosenbrock23())

# Access results
println("Gas phase NH3: ", sol[sys.g.NH3.M][end], " mol/m³")
println("Aerosol water: ", sol[sys.aq.W][end], " kg/m³")
```

## Analysis

This implementation provides a comprehensive thermodynamic equilibrium model suitable for atmospheric chemistry applications, with full validation against the original FORTRAN reference code.
