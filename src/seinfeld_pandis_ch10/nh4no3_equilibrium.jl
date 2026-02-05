"""
Ammonium nitrate equilibrium equations from Seinfeld & Pandis (2006) Chapter 10.

NH4NO3 is a semi-volatile aerosol species whose gas-particle partitioning
depends strongly on temperature and relative humidity.
"""

"""
    NH4NO3Equilibrium(; name=:NH4NO3Equilibrium)

A component that calculates NH4NO3 gas-aerosol equilibrium.

Implements the equilibrium reactions and constants from Section 10.4.3 of
Seinfeld & Pandis (2006).

## Solid NH4NO3 equilibrium (below DRH)

Eq. 10.87: NH3(g) + HNO3(g) ⇌ NH4NO3(s)

Eq. 10.91 - Dissociation constant (ppb² units):

```math
\\ln K_p = 84.6 - \\frac{24220}{T} - 6.1 \\ln\\frac{T}{298}
```

## Aqueous NH4NO3 equilibrium (above DRH)

Eq. 10.92: NH3(g) + HNO3(g) ⇌ NH4+ + NO3-

Eq. 10.97 - Equilibrium constant:

```math
K_{AN} = 4 \\times 10^{17} \\exp\\left\\{64.7\\left(\\frac{298}{T} - 1\\right) + 11.51\\left[1 + \\ln\\frac{298}{T} - \\frac{298}{T}\\right]\\right\\}
```

(units: mol² kg⁻² atm⁻²)

## Ionic strength fraction (Eq. 10.100)

```math
Y = \\frac{[NH_4NO_3]}{[NH_4NO_3] + 3[(NH_4)_2SO_4]}
```

# Parameters

  - `T`: Temperature (K)
  - `RH`: Relative humidity (0-1)
  - `C_NH42SO4`: (NH4)2SO4 concentration for ionic strength calculation (mol/m³)

# Variables

  - `ln_Kp`: Log of solid dissociation constant (dimensionless for ppb² units)
  - `Kp`: Solid dissociation constant (ppb²)
  - `K_AN`: Aqueous equilibrium constant (mol² kg⁻² atm⁻²)
  - `DRH`: Deliquescence RH of NH4NO3 (dimensionless)
  - `is_aqueous`: 1 if RH > DRH (aqueous phase), 0 otherwise (dimensionless)
  - `Y`: Ionic strength fraction (dimensionless)

# Example

        # Temperature coefficients with units for Eq. 10.91 and 10.88

```julia
using ModelingToolkit, Aerosol
@mtkbuild sys = NH4NO3Equilibrium()
```
"""
@component function NH4NO3Equilibrium(; name = :NH4NO3Equilibrium)
    @constants begin
        T_ref = 298.0, [description = "Reference temperature", unit = u"K"]
        # Temperature coefficients with units for Eq. 10.91 and 10.88
        c1_Kp = 84.6,
        [description = "First coefficient in Eq. 10.91 (dimensionless)", unit = u"1"]
        c2_Kp = 24220.0, [description = "Second coefficient in Eq. 10.91", unit = u"K"]
        c3_Kp = 6.1,
        [description = "Third coefficient in Eq. 10.91 (dimensionless)", unit = u"1"]
        c1_DRH = 723.7, [description = "DRH coefficient", unit = u"K"]
        c2_DRH = 1.6954, [description = "DRH constant (dimensionless)", unit = u"1"]
    end

    @parameters begin
        T = 298.0, [description = "Temperature", unit = u"K"]
        RH = 0.5, [description = "Relative humidity (0-1) (dimensionless)", unit = u"1"]
    end

    @variables begin
        ln_Kp(t),
        [description = "Log of solid dissociation constant (ppb² basis) (dimensionless)",
            unit = u"1"]
        DRH_out(t), [description = "Deliquescence RH (dimensionless)", unit = u"1"]
        is_aqueous(t),
        [description = "Phase indicator (1=aqueous, 0=solid) (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 10.91: Solid NH4NO3 dissociation constant (ln(Kp) where Kp is in ppb²)
        ln_Kp ~ c1_Kp - c2_Kp/T - c3_Kp * log(T/T_ref),  # Eq. 10.91

        # Eq. 10.88: NH4NO3 DRH vs temperature
        DRH_out ~ exp(c1_DRH/T + c2_DRH) / 100.0,  # Eq. 10.88

        # Phase determination: aqueous if RH > DRH
        is_aqueous ~ ifelse(RH > DRH_out, 1.0, 0.0)
    ]

    return System(eqs, t; name)
end

"""
    nh4no3_Kp(T)

Calculate the solid NH4NO3 dissociation constant at temperature T.

Implements Equation 10.91 from Seinfeld & Pandis (2006):

```math
\\ln K_p = 84.6 - \\frac{24220}{T} - 6.1 \\ln\\frac{T}{298}
```

# Arguments

  - `T`: Temperature (K)

# Returns

Kp in ppb² units
"""
function nh4no3_Kp(T::Real)
    return exp(84.6 - 24220/T - 6.1 * log(T/298.0))
end

"""
    nh4no3_K_AN(T)

Calculate the aqueous NH4NO3 equilibrium constant at temperature T.

Implements Equation 10.97 from Seinfeld & Pandis (2006).

# Arguments

  - `T`: Temperature (K)

# Returns

K_AN in mol² kg⁻² atm⁻² units
"""
function nh4no3_K_AN(T::Real)
    T_ref = 298.0
    K_AN_ref = 4.0e17  # mol² kg⁻² atm⁻²
    a = 64.7
    b = 11.51

    return K_AN_ref * exp(a * (T_ref/T - 1) + b * (1 + log(T_ref/T) - T_ref/T))
end

"""
    ionic_strength_fraction(C_NH4NO3, C_NH42SO4)

Calculate the ionic strength fraction Y (Eq. 10.100).

# Arguments

  - `C_NH4NO3`: NH4NO3 concentration (any units, same as C_NH42SO4)
  - `C_NH42SO4`: (NH4)2SO4 concentration (any units, same as C_NH4NO3)

# Returns

Ionic strength fraction Y (dimensionless, 0-1)
"""
function ionic_strength_fraction(C_NH4NO3::Real, C_NH42SO4::Real)
    return C_NH4NO3 / (C_NH4NO3 + 3 * C_NH42SO4 + 1e-20)
end

"""
    equilibrium_constant_temperature(K_ref, T, T_ref, a, b)

Calculate temperature-dependent equilibrium constant using the standard form.

```math
K(T) = K(298) \\exp\\left\\{a\\left(\\frac{298}{T} - 1\\right) + b\\left[1 + \\ln\\frac{298}{T} - \\frac{298}{T}\\right]\\right\\}
```

This form is used for equilibrium constants in Table 10.7 of Seinfeld & Pandis (2006).

# Arguments

  - `K_ref`: Equilibrium constant at reference temperature
  - `T`: Temperature (K)
  - `T_ref`: Reference temperature (K), typically 298
  - `a`: First temperature coefficient
  - `b`: Second temperature coefficient

# Returns

Equilibrium constant at temperature T
"""
function equilibrium_constant_temperature(
        K_ref::Real, T::Real, T_ref::Real, a::Real, b::Real)
    return K_ref * exp(a * (T_ref/T - 1) + b * (1 + log(T_ref/T) - T_ref/T))
end

# Table 10.7: Equilibrium constant parameters
# Format: (K_298, a, b, units_description)
const EQUILIBRIUM_CONSTANTS = Dict{Symbol, NTuple{4, Any}}(
    :NaCl_HNO3 => (3.96, 5.50, -2.18, "dimensionless"),
    :HSO4_dissoc => (1.01e-2, 8.85, 25.14, "mol/kg"),
    :NH3_HNO3_aq => (4.0e17, 64.7, 11.51, "mol² kg⁻² atm⁻²"),
    :HCl_dissoc => (2.03e6, 30.21, 19.91, "mol² kg⁻² atm⁻¹"),
    :NH3_HCl_aq => (2.12e17, 65.08, 14.51, "mol² kg⁻² atm⁻²"),
    :Na2SO4_dissoc => (0.48, 0.98, 39.57, "mol³ kg⁻³"),
    :NH42SO4_dissoc => (1.425, -2.65, 38.55, "mol³ kg⁻³"),
    :HNO3_dissoc => (3.638e6, 29.47, 16.84, "mol² kg⁻² atm⁻¹"),
    :NH4Cl_dissoc => (1.039e-16, -71.04, 2.40, "atm²"),
    :NH4NO3_solid => (3.35e16, 75.11, -13.5, "atm⁻²"),
    :NaCl_dissoc => (37.74, -1.57, 16.89, "mol² kg⁻²"),
    :NaHSO4_dissoc => (2.44e4, 0.79, 4.53, "mol² kg⁻²"),
    :NaNO3_dissoc => (11.97, -8.22, 16.0, "mol² kg⁻²")
)

"""
    get_equilibrium_constant(reaction::Symbol, T::Real)

Get the equilibrium constant for a reaction at temperature T using Table 10.7 data.

# Arguments

  - `reaction`: Symbol identifying the reaction (e.g., :NH4NO3_solid, :HSO4_dissoc)
  - `T`: Temperature (K)

# Returns

Equilibrium constant at temperature T

# Available reactions

  - :NaCl_HNO3 - NaCl(s) + HNO3(g) ⇌ NaNO3(s) + HCl(g)
  - :HSO4_dissoc - HSO4⁻ ⇌ H⁺ + SO4²⁻
  - :NH3_HNO3_aq - NH3(g) + HNO3(g) ⇌ NH4⁺ + NO3⁻
  - :HCl_dissoc - HCl(g) ⇌ H⁺ + Cl⁻
  - :NH3_HCl_aq - NH3(g) + HCl(g) ⇌ NH4⁺ + Cl⁻
  - :Na2SO4_dissoc - Na2SO4(s) ⇌ 2Na⁺ + SO4²⁻
  - :NH42SO4_dissoc - (NH4)2SO4(s) ⇌ 2NH4⁺ + SO4²⁻
  - :HNO3_dissoc - HNO3(g) ⇌ H⁺ + NO3⁻
  - :NH4Cl_dissoc - NH4Cl(s) ⇌ NH3(g) + HCl(g)
  - :NH4NO3_solid - NH3(g) + HNO3(g) ⇌ NH4NO3(s)
  - :NaCl_dissoc - NaCl(s) ⇌ Na⁺ + Cl⁻
  - :NaHSO4_dissoc - NaHSO4(s) ⇌ Na⁺ + HSO4⁻
  - :NaNO3_dissoc - NaNO3(s) ⇌ Na⁺ + NO3⁻
"""
function get_equilibrium_constant(reaction::Symbol, T::Real)
    if !haskey(EQUILIBRIUM_CONSTANTS, reaction)
        error("Unknown reaction: $reaction. Available: $(keys(EQUILIBRIUM_CONSTANTS))")
    end

    K_ref, a, b, _ = EQUILIBRIUM_CONSTANTS[reaction]
    return equilibrium_constant_temperature(K_ref, T, 298.0, a, b)
end
