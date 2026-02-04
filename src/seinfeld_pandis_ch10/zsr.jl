"""
ZSR (Zdanovskii-Stokes-Robinson) water content equation from Seinfeld & Pandis (2006).

The ZSR mixing rule calculates the aerosol liquid water content for multicomponent
solutions based on the water uptake of individual salts.
"""

"""
    ZSRWaterContent(; name=:ZSRWaterContent, n_species=3)

A component that calculates aerosol liquid water content using the ZSR mixing rule.

Implements Equation 10.98 from Seinfeld & Pandis (2006):

```math
W = \\sum_i \\frac{C_i}{m_{i,0}(a_w)}
```

where:

  - `W` is the total aerosol liquid water content (kg water per m³ air)
  - `Cᵢ` is the concentration of species i (mol/m³ air)
  - `mᵢ,₀(aᵤ)` is the molality of species i in a binary solution at water activity aᵤ

The water activity equals the relative humidity for equilibrium aerosols (Eq. 10.63).

# Parameters

  - `RH`: Relative humidity as a fraction (0-1)
  - `C[1:n_species]`: Concentrations of each salt species (mol/m³)
  - `m0[1:n_species]`: Molalities of each species at saturation in binary solution (mol/kg)

# Variables

  - `W`: Total aerosol liquid water content (kg/m³)
  - `α_w`: Water activity (dimensionless)

# Example

```julia
using ModelingToolkit, Aerosol
@mtkbuild sys = ZSRWaterContent(n_species = 2)
```
"""
@component function ZSRWaterContent(; name = :ZSRWaterContent, n_species::Int = 3)
    @parameters begin
        RH = 0.8, [description = "Relative humidity (0-1) (dimensionless)", unit = u"1"]
        C[1:n_species] = fill(1e-6, n_species),
        [description = "Salt concentrations", unit = u"mol/m^3"]
        m0[1:n_species] = fill(5.0, n_species),
        [description = "Binary solution molalities at water activity α_w", unit = u"mol/kg"]
    end

    @variables begin
        W(t), [description = "Aerosol liquid water content", unit = u"kg/m^3"]
        α_w(t), [description = "Water activity (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 10.63: Water activity equals RH for equilibrium aerosol
        α_w ~ RH,  # Eq. 10.63

        # Eq. 10.98: ZSR mixing rule for water content
        W ~ sum(C[i] / m0[i] for i in 1:n_species)  # Eq. 10.98
    ]

    return System(eqs, t; name)
end

# Binary solution molality data at various water activities
# Data from Table 10.2 and empirical fits in Section 10.2

"""
    binary_molality_nh42so4(a_w)

Calculate the molality of (NH4)2SO4 in binary aqueous solution at water activity a_w.

Based on polynomial fits to water activity data.
"""
function binary_molality_nh42so4(a_w)
    # Empirical fit for (NH4)2SO4 binary solution
    # Valid for a_w between ~0.4 and 1.0
    if a_w >= 0.99
        return 0.01
    elseif a_w <= 0.40
        return 29.0  # Near saturation
    end
    # Polynomial fit: m = f(a_w)
    return 5.92 * (1 - a_w)^0.45 + 0.1
end

"""
    binary_molality_nh4no3(a_w)

Calculate the molality of NH4NO3 in binary aqueous solution at water activity a_w.
"""
function binary_molality_nh4no3(a_w)
    if a_w >= 0.99
        return 0.01
    elseif a_w <= 0.30
        return 30.0  # Near saturation
    end
    return 12.0 * (1 - a_w)^0.52 + 0.1
end

"""
    binary_molality_nacl(a_w)

Calculate the molality of NaCl in binary aqueous solution at water activity a_w.
"""
function binary_molality_nacl(a_w)
    if a_w >= 0.99
        return 0.01
    elseif a_w <= 0.45
        return 14.0  # Near saturation
    end
    return 6.0 * (1 - a_w)^0.43 + 0.05
end

"""
    zsr_water_content(RH, concentrations::Dict{Symbol,Float64})

Calculate aerosol liquid water content using the ZSR mixing rule (Eq. 10.98).

# Arguments

  - `RH`: Relative humidity as a fraction (0-1)
  - `concentrations`: Dictionary mapping salt symbols to concentrations (mol/m³)

# Returns

Aerosol liquid water content W (kg/m³)

# Supported salts

  - :NH42SO4 - Ammonium sulfate
  - :NH4NO3 - Ammonium nitrate
  - :NaCl - Sodium chloride

# Example

```julia
W = zsr_water_content(0.8, Dict(:NH42SO4 => 1e-6, :NH4NO3 => 2e-6))  # Eq. 10.63
```
"""
function zsr_water_content(RH::Real, concentrations::Dict{Symbol, Float64})
    a_w = RH  # Eq. 10.63

    # Map of salts to their binary molality functions
    molality_funcs = Dict{Symbol, Function}(
        :NH42SO4 => binary_molality_nh42so4,
        :NH4NO3 => binary_molality_nh4no3,
        :NaCl => binary_molality_nacl
    )

    W = 0.0
    for (salt, C) in concentrations
        if haskey(molality_funcs, salt)
            m0 = molality_funcs[salt](a_w)
            W += C / m0
        end
    end

    return W
end
