"""
Kelvin effect equations from Seinfeld & Pandis (2006) Chapter 10.

The Kelvin effect describes the enhancement of vapor pressure over a curved
liquid surface relative to a flat surface. This is particularly important
for small droplets where surface curvature is significant.
"""

"""
    KelvinEffect(; name=:KelvinEffect)

A component that calculates the Kelvin effect on vapor pressure over curved droplet surfaces.

Implements Equation 10.86 from Seinfeld & Pandis (2006):

```math
p_A = p°_A \\exp\\left(\\frac{2\\sigma M}{RT \\rho_l R_p}\\right)
```

where:
- `p_A` is the vapor pressure over the droplet
- `p°_A` is the saturation vapor pressure over a flat surface
- `σ` is the surface tension
- `M` is the molar mass
- `R` is the gas constant
- `T` is temperature
- `ρ_l` is the liquid density
- `R_p` is the droplet radius

# Variables
- `S`: Saturation ratio (p_A / p°_A) (dimensionless)
- `ln_S`: Natural log of saturation ratio (dimensionless)

# Parameters
- `T`: Temperature (K)
- `R_p`: Particle radius (m)
- `σ`: Surface tension (N/m)
- `M_mol`: Molar mass (kg/mol)
- `ρ_l`: Liquid density (kg/m³)

# Example
```julia
using ModelingToolkit, Aerosol
@mtkbuild sys = KelvinEffect()
```
"""
@component function KelvinEffect(; name=:KelvinEffect)
    @constants begin
        R_gas = 8.314, [description = "Universal gas constant", unit = u"J/mol/K"]
    end

    @parameters begin
        T = 298.15, [description = "Temperature", unit = u"K"]
        R_p = 50e-9, [description = "Particle radius", unit = u"m"]
        σ = 0.072, [description = "Surface tension", unit = u"N/m"]
        M_mol = 0.018, [description = "Molar mass", unit = u"kg/mol"]
        ρ_l = 1000.0, [description = "Liquid density", unit = u"kg/m^3"]
    end

    @variables begin
        ln_S(t), [description = "Natural log of saturation ratio (dimensionless)", unit = u"1"]
        S(t), [description = "Saturation ratio p_A/p°_A (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 10.86: Kelvin equation
        ln_S ~ 2 * σ * M_mol / (R_gas * T * ρ_l * R_p),  # Eq. 10.86
        S ~ exp(ln_S),  # Saturation ratio
    ]

    return System(eqs, t; name)
end


"""
    kelvin_saturation_ratio(T, R_p, σ, M_mol, ρ_l)

Calculate the Kelvin saturation ratio for a droplet.

Implements Equation 10.86 from Seinfeld & Pandis (2006).

# Arguments
- `T`: Temperature (K)
- `R_p`: Particle radius (m)
- `σ`: Surface tension (N/m)
- `M_mol`: Molar mass (kg/mol)
- `ρ_l`: Liquid density (kg/m³)

# Returns
Saturation ratio S = p_A/p°_A (dimensionless)

# Example
For a 50 nm water droplet at 298 K:
```julia
S = kelvin_saturation_ratio(298.0, 50e-9, 0.072, 0.018, 1000.0)
# S ≈ 1.024 (2.4% vapor pressure enhancement)
```
"""
function kelvin_saturation_ratio(T, R_p, σ, M_mol, ρ_l)
    R_gas = 8.314  # J/(mol·K)
    return exp(2 * σ * M_mol / (R_gas * T * ρ_l * R_p))
end
