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

Note: Units are documented in descriptions but not enforced via DynamicQuantities
due to non-SI units (atm) used in atmospheric chemistry conventions.
All concentrations are in mol/L (M), pressures in atm, temperatures in K.
"""

# =============================================================================
# Physical Constants
# =============================================================================

"""
Gas constant in atm L mol^-1 K^-1 for Henry's law calculations
"""
const R_GAS_ATM = 0.08205  # atm L mol^-1 K^-1

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
# Units: M atm^-1
# =============================================================================

const HENRY_CONSTANTS_298 = Dict(
    :O2 => 1.3e-3,
    :O3 => 1.1e-2,
    :CO2 => 3.4e-2,
    :SO2 => 1.23,
    :NH3 => 62.0,
    :H2O2 => 1.0e5,
    :HNO3 => 2.1e5,
    :HCHO => 2.5,  # Without diol formation
    :HCHO_diol => 6.3e3,  # With diol formation
    :HCOOH => 3.6e3,
    :CH3OOH => 310.0,
)

# =============================================================================
# Heat of Dissolution (Table 7.3)
# Units: kcal mol^-1
# =============================================================================

const DELTA_H_DISSOLUTION = Dict(
    :CO2 => -4.85,
    :NH3 => -8.17,
    :SO2 => -6.25,
    :H2O2 => -14.5,
    :O3 => -5.04,
)

# =============================================================================
# Basic Henry's Law Component (Fixed Temperature)
# =============================================================================

"""
    HenrysLaw(; name=:HenrysLaw)

Basic Henry's law equilibrium at fixed temperature.

Implements Eq 7.3: [A(aq)] = H_A * p_A

Parameters:
- H_A: Henry's law constant (M atm^-1)

Variables:
- p_A: Partial pressure of gas A (atm)
- C_aq: Aqueous concentration (mol/L = M)

This component establishes the equilibrium relationship between gas-phase
partial pressure and aqueous-phase concentration.
"""
@component function HenrysLaw(; name=:HenrysLaw)
    @parameters begin
        H_A, [description = "Henry's law constant at reference temperature (M/atm)"]
    end

    @variables begin
        p_A(t), [description = "Partial pressure of gas species A (atm)"]
        C_aq(t), [description = "Aqueous concentration of species A (mol/L)"]
    end

    eqs = [
        # Eq 7.3: Henry's law equilibrium
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
- H_298: Henry's law constant at 298 K (M atm^-1)
- dH_diss: Heat of dissolution (J mol^-1)

Variables:
- T: Temperature (K)
- p_A: Partial pressure (atm)
- C_aq: Aqueous concentration (M)
- H_T: Temperature-corrected Henry's law constant (M atm^-1)

Also calculates:
- w_L: Liquid water mixing ratio (vol/vol)
- f_A: Distribution factor (dimensionless)
- X_aq: Aqueous fraction (dimensionless)
"""
@component function HenrysLawTemperature(; name=:HenrysLawTemp)
    @constants begin
        R_gas = 8.314, [description = "Gas constant (J/mol/K)"]
        T_ref = 298.0, [description = "Reference temperature (K)"]
        R_atm = 0.08205, [description = "Gas constant for distribution factor (L*atm/mol/K)"]
    end

    @parameters begin
        H_298, [description = "Henry's law constant at 298 K (M/atm)"]
        dH_diss, [description = "Heat of dissolution (J/mol)"]
        L, [description = "Liquid water content (g/m^3)"]
    end

    @variables begin
        T(t), [description = "Temperature (K)"]
        p_A(t), [description = "Partial pressure of gas species A (atm)"]
        C_aq(t), [description = "Aqueous concentration of species A (mol/L)"]
        H_T(t), [description = "Temperature-corrected Henry's law constant (M/atm)"]
        w_L(t), [description = "Liquid water mixing ratio (vol/vol, dimensionless)"]
        f_A(t), [description = "Distribution factor (dimensionless)"]
        X_aq(t), [description = "Aqueous fraction (dimensionless)"]
    end

    eqs = [
        # Eq 7.5: Temperature dependence of Henry's law constant
        # H_A(T) = H_A(T_ref) * exp((dH/R) * (1/T_ref - 1/T))
        H_T ~ H_298 * exp((dH_diss / R_gas) * (1/T_ref - 1/T)),

        # Eq 7.3: Henry's law equilibrium at temperature T
        C_aq ~ H_T * p_A,

        # Eq 7.1: Liquid water mixing ratio
        # w_L (vol water/vol air) = 10^-6 * L (g m^-3)
        # Assuming density of water = 10^6 g/m^3
        w_L ~ 1e-6 * L,

        # Eq 7.7: Distribution factor
        # f_A = H_A * R * T * w_L
        f_A ~ H_T * R_atm * T * w_L,

        # Eq 7.9: Aqueous fraction
        # X_aq = f_A / (1 + f_A)
        X_aq ~ f_A / (1 + f_A),
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
- H_intrinsic: Intrinsic Henry's law constant (M atm^-1)
- K_1: First dissociation constant (M)
- K_2: Second dissociation constant (M), set to 0 for monoprotic acids

Variables:
- H_plus: Hydrogen ion concentration (M)
- H_eff: Effective Henry's law constant (M atm^-1)
"""
@component function EffectiveHenrysLaw(; name=:EffectiveHenrys)
    @parameters begin
        H_intrinsic, [description = "Intrinsic Henry's law constant (M/atm)"]
        K_1, [description = "First dissociation constant (M)"]
        K_2, [description = "Second dissociation constant (M)"]
    end

    @variables begin
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
        H_eff(t), [description = "Effective Henry's law constant (M/atm)"]
    end

    eqs = [
        # General form for diprotic species
        # H* = H * (1 + K_1/[H+] + K_1*K_2/[H+]^2)
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
- H: Intrinsic Henry's law constant (M atm^-1)
- K1: First dissociation constant (M)
- K2: Second dissociation constant (M), use 0 for monoprotic
- H_plus: Hydrogen ion concentration (M)

Returns:
- H_eff: Effective Henry's law constant (M atm^-1)
"""
function effective_henrys_constant(H, K1, K2, H_plus)
    return H * (1 + K1/H_plus + K1*K2/H_plus^2)
end

"""
    aqueous_fraction(H, R, T, w_L)

Calculate the fraction of species in aqueous phase (Eq 7.9).

Arguments:
- H: Henry's law constant (M atm^-1)
- R: Gas constant (0.08205 atm L mol^-1 K^-1)
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
- H: Henry's law constant (M atm^-1)
- R: Gas constant (0.08205 atm L mol^-1 K^-1)
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
- H_298: Henry's law constant at 298 K (M atm^-1)
- dH_diss: Heat of dissolution (J mol^-1)
- T: Temperature (K)

Returns:
- H_T: Henry's law constant at temperature T (M atm^-1)
"""
function henrys_constant_at_T(H_298, dH_diss, T)
    R = 8.314  # J mol^-1 K^-1
    T_ref = 298.0  # K
    return H_298 * exp((dH_diss / R) * (1/T_ref - 1/T))
end
