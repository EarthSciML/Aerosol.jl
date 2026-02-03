"""
    Henry's Law Components

Implements gas-liquid equilibrium calculations from Seinfeld & Pandis Chapter 7,
Section 7.2.

Key equations:
- Eq 7.1: w_L = 10^-6 * L (liquid water mixing ratio)
- Eq 7.3: [A(aq)] = H_A * p_A (Henry's law)
- Eq 7.5: H_A(T2) = H_A(T1) * exp((ΔH_A/R) * (1/T1 - 1/T2)) (temperature dependence)
- Eq 7.7: f_A = H_A * R * T * w_L (distribution factor)
- Eq 7.9: X_aq = (H_A * R * T * w_L) / (1 + H_A * R * T * w_L) (aqueous fraction)

All pressures are in Pa (SI). Henry's law constants from the textbook (M/atm) are
converted to mol/m³/Pa for use in @component functions (SI units).
Module-level constants retain M/Pa (mol/L/Pa) units for use in utility functions.
Concentrations in @component functions are in mol/m³ (SI), temperatures in K.
"""

# =============================================================================
# Physical Constants
# =============================================================================

"""
Standard atmosphere in Pa, used for converting Henry's law constants from M/atm to M/Pa.
"""
const ATM_TO_PA = 101325.0  # Pa/atm

"""
Gas constant in Pa L mol^-1 K^-1 for Henry's law calculations.
R = 8.314 J/(mol·K) = 8.314 Pa·m³/(mol·K) = 8314.0 Pa·L/(mol·K) / 1000
Actually R in L·atm/(mol·K) = 0.08205; in L·Pa/(mol·K) = 0.08205 * 101325 = 8314.46
"""
const R_GAS_PA = 8314.46  # Pa L mol^-1 K^-1

"""
Gas constant in J mol^-1 K^-1 for van't Hoff equation
"""
const R_GAS_J = 8.314  # J mol^-1 K^-1

"""
Reference temperature for Henry's law constants (K)
"""
const T_REF = 298.0  # K

# =============================================================================
# Henry's Law Constants at 298 K (Table 7.2)
# Units: M/Pa (converted from M/atm by dividing by 101325)
# =============================================================================

const HENRY_CONSTANTS_298 = Dict(
    :O2 => 1.3e-3 / ATM_TO_PA,
    :O3 => 1.1e-2 / ATM_TO_PA,
    :CO2 => 3.4e-2 / ATM_TO_PA,
    :SO2 => 1.23 / ATM_TO_PA,
    :NH3 => 62.0 / ATM_TO_PA,
    :H2O2 => 1.0e5 / ATM_TO_PA,
    :HNO3 => 2.1e5 / ATM_TO_PA,
    :HCHO => 2.5 / ATM_TO_PA,  # Without diol formation
    :HCHO_diol => 6.3e3 / ATM_TO_PA,  # With diol formation
    :HCOOH => 3.6e3 / ATM_TO_PA,
    :CH3OOH => 310.0 / ATM_TO_PA
)

# =============================================================================
# Heat of Dissolution (Table 7.3)
# Units: J mol^-1 (converted from kcal mol^-1 by multiplying by 4184)
# =============================================================================

const DELTA_H_DISSOLUTION = Dict(
    :CO2 => -4.85 * 4184,
    :NH3 => -8.17 * 4184,
    :SO2 => -6.25 * 4184,
    :H2O2 => -14.5 * 4184,
    :O3 => -5.04 * 4184
)

# =============================================================================
# Basic Henry's Law Component (Fixed Temperature)
# =============================================================================

"""
    HenrysLaw(; name=:HenrysLaw)

Basic Henry's law equilibrium at fixed temperature.

Implements Eq 7.3: [A(aq)] = H_A * p_A

Parameters:

  - H_A: Henry's law constant (mol/m³/Pa)

Variables:

  - p_A: Partial pressure of gas A (Pa)
  - C_aq: Aqueous concentration (mol/m³)

This component establishes the equilibrium relationship between gas-phase
partial pressure and aqueous-phase concentration.
"""
@component function HenrysLaw(; name = :HenrysLaw)
    @parameters begin
        H_A,
        [
            description = "Henry's law constant at reference temperature", unit = u"mol/m^3/Pa"]
    end

    @variables begin
        p_A(t), [description = "Partial pressure of gas species A", unit = u"Pa"]
        C_aq(t), [description = "Aqueous concentration of species A", unit = u"mol/m^3"]
    end

    eqs = [
    # Eq 7.3: Henry's law equilibrium
    # (mol/m^3/Pa) * (Pa) = (mol/m^3) ✓
        C_aq ~ H_A * p_A,
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Temperature-Dependent Henry's Law Component
# =============================================================================

"""
    HenrysLawTemperature(; name=:HenrysLawTemp)

Temperature-dependent Henry's law equilibrium.

Implements:

  - Eq 7.3: [A(aq)] = H_A(T) * p_A
  - Eq 7.5: H_A(T) = H_A(T_ref) * exp((ΔH_A/R) * (1/T_ref - 1/T))

Parameters:

  - H_298: Henry's law constant at 298 K (mol/m³/Pa)
  - dH_diss: Heat of dissolution (J/mol)

Variables:

  - T: Temperature (K)
  - p_A: Partial pressure (Pa)
  - C_aq: Aqueous concentration (mol/m³)
  - H_T: Temperature-corrected Henry's law constant (mol/m³/Pa)

Also calculates:

  - w_L: Liquid water mixing ratio (dimensionless, vol/vol)
  - f_A: Distribution factor (dimensionless)
  - X_aq: Aqueous fraction (dimensionless)
"""
@component function HenrysLawTemperature(; name = :HenrysLawTemp)
    @constants begin
        R_gas = 8.314, [description = "Gas constant", unit = u"J/mol/K"]
        T_ref = 298.0, [description = "Reference temperature", unit = u"K"]
        # Inverse water density: 1/(10^6 g/m^3) = 10^-6 m^3/g
        rho_water_inv = 1e-6,
        [description = "Inverse water density for LWC conversion", unit = u"m^3/g"]
    end

    @parameters begin
        H_298, [description = "Henry's law constant at 298 K", unit = u"mol/m^3/Pa"]
        dH_diss, [description = "Heat of dissolution", unit = u"J/mol"]
        L, [description = "Liquid water content", unit = u"g/m^3"]
    end

    @variables begin
        T(t), [description = "Temperature", unit = u"K"]
        p_A(t), [description = "Partial pressure of gas species A", unit = u"Pa"]
        C_aq(t), [description = "Aqueous concentration of species A", unit = u"mol/m^3"]
        H_T(t),
        [description = "Temperature-corrected Henry's law constant", unit = u"mol/m^3/Pa"]
        w_L(t), [description = "Liquid water mixing ratio (dimensionless)", unit = u"1"]
        f_A(t), [description = "Distribution factor (dimensionless)", unit = u"1"]
        X_aq(t), [description = "Aqueous fraction (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq 7.5: Temperature dependence of Henry's law constant
        # H_A(T) = H_A(T_ref) * exp((dH/R) * (1/T_ref - 1/T))
        # Units: (mol/m^3/Pa) * exp((J/mol) / (J/mol/K) * (1/K)) = (mol/m^3/Pa) ✓
        H_T ~ H_298 * exp((dH_diss / R_gas) * (1/T_ref - 1/T)),

        # Eq 7.3: Henry's law equilibrium at temperature T
        # Units: (mol/m^3/Pa) * (Pa) = (mol/m^3) ✓
        C_aq ~ H_T * p_A,

        # Eq 7.1: Liquid water mixing ratio
        # w_L (vol water/vol air) = L / rho_water = L * rho_water_inv
        # Units: (g/m^3) * (m^3/g) = (1) ✓
        w_L ~ rho_water_inv * L,

        # Eq 7.7: Distribution factor
        # f_A = H_A * R * T * w_L
        # Using R_gas in J/(mol·K) = Pa·m³/(mol·K), consistent with H_T in mol/m³/Pa
        # Units: (mol/m^3/Pa) * (Pa*m^3/mol/K) * (K) * (1) = (1) ✓
        f_A ~ H_T * R_gas * T * w_L,

        # Eq 7.9: Aqueous fraction
        # X_aq = f_A / (1 + f_A)
        X_aq ~ f_A / (1 + f_A)
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Effective Henry's Law Constants (with dissociation)
# =============================================================================

"""
    EffectiveHenrysLaw(; name=:EffectiveHenrys)

Effective Henry's law constant accounting for aqueous dissociation.

For species that dissociate in water, the effective Henry's law constant H*
is larger than the intrinsic Henry's law constant H because dissociation
products contribute to the total dissolved amount.

Implements:

  - Eq 7.24: H*_CO2 = H_CO2 * (1 + K_c1/[H+] + K_c1*K_c2/[H+]^2)
  - Eq 7.41: H*_SO2 = H_SO2 * (1 + K_s1/[H+] + K_s1*K_s2/[H+]^2)
  - Eq 7.61: H*_HNO3 ~ H_HNO3 * K_n1 / [H+]

Parameters:

  - H_intrinsic: Intrinsic Henry's law constant (mol/m³/Pa)
  - K_1: First dissociation constant (mol/m³)
  - K_2: Second dissociation constant (mol/m³), set to 0 for monoprotic acids

Variables:

  - H_plus: Hydrogen ion concentration (mol/m³)
  - H_eff: Effective Henry's law constant (mol/m³/Pa)
"""
@component function EffectiveHenrysLaw(; name = :EffectiveHenrys)
    @parameters begin
        H_intrinsic, [description = "Intrinsic Henry's law constant", unit = u"mol/m^3/Pa"]
        K_1, [description = "First dissociation constant", unit = u"mol/m^3"]
        K_2, [description = "Second dissociation constant", unit = u"mol/m^3"]
    end

    @variables begin
        H_plus(t), [description = "Hydrogen ion concentration", unit = u"mol/m^3"]
        H_eff(t), [description = "Effective Henry's law constant", unit = u"mol/m^3/Pa"]
    end

    eqs = [
    # General form for diprotic species
    # H* = H * (1 + K_1/[H+] + K_1*K_2/[H+]^2)
    # Units: K_1/H_plus = (mol/m^3)/(mol/m^3) = (1) ✓
    # Units: K_1*K_2/H_plus^2 = (mol/m^3)^2/(mol/m^3)^2 = (1) ✓
        H_eff ~ H_intrinsic * (1 + K_1/H_plus + K_1*K_2/H_plus^2),
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Utility Functions
# =============================================================================

"""
    effective_henrys_constant(H, K1, K2, H_plus)

Calculate effective Henry's law constant accounting for dissociation.

Arguments:

  - H: Intrinsic Henry's law constant (mol/L/Pa)
  - K1: First dissociation constant (mol/L)
  - K2: Second dissociation constant (mol/L), use 0 for monoprotic
  - H_plus: Hydrogen ion concentration (mol/L)

Returns:

  - H_eff: Effective Henry's law constant (mol/L/Pa)
"""
function effective_henrys_constant(H, K1, K2, H_plus)
    return H * (1 + K1/H_plus + K1*K2/H_plus^2)
end

"""
    aqueous_fraction(H, R, T, w_L)

Calculate the fraction of species in aqueous phase (Eq 7.9).

Arguments:

  - H: Henry's law constant (mol/L/Pa)
  - R: Gas constant (8314.46 Pa L mol^-1 K^-1)
  - T: Temperature (K)
  - w_L: Liquid water mixing ratio (vol/vol)

Returns:

  - X_aq: Fraction in aqueous phase (dimensionless)
"""
function aqueous_fraction(H, R, T, w_L)
    f = H * R * T * w_L
    return f / (1 + f)
end

"""
    distribution_factor(H, R, T, w_L)

Calculate the distribution factor f_A (Eq 7.7).

Arguments:

  - H: Henry's law constant (mol/L/Pa)
  - R: Gas constant (8314.46 Pa L mol^-1 K^-1)
  - T: Temperature (K)
  - w_L: Liquid water mixing ratio (vol/vol)

Returns:

  - f_A: Distribution factor (dimensionless)
"""
function distribution_factor(H, R, T, w_L)
    return H * R * T * w_L
end

"""
    henrys_constant_at_T(H_298, dH_diss, T)

Calculate Henry's law constant at temperature T using van't Hoff equation (Eq 7.5).

Arguments:

  - H_298: Henry's law constant at 298 K (mol/L/Pa)
  - dH_diss: Heat of dissolution (J/mol)
  - T: Temperature (K)

Returns:

  - H_T: Henry's law constant at temperature T (mol/L/Pa)
"""
function henrys_constant_at_T(H_298, dH_diss, T)
    R = 8.314  # J mol^-1 K^-1
    T_ref = 298.0  # K
    return H_298 * exp((dH_diss / R) * (1/T_ref - 1/T))
end
