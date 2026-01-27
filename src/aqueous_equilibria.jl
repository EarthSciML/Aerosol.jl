"""
    Aqueous Equilibria Components

Implements dissociation equilibria for atmospheric species in cloud droplets
from Seinfeld & Pandis Chapter 7, Section 7.3.

Includes:
- Water autoionization (Eq 7.10-7.13)
- CO2 system (Section 7.3.2, Eq 7.14-7.24)
- SO2 system (Section 7.3.3, Eq 7.30-7.45)
- NH3 system (Section 7.3.4, Eq 7.46-7.51)
- HNO3 system (Section 7.3.5, Eq 7.53-7.61)
- H2O2 system (Section 7.3.6)

Temperature dependence follows Table 7.4 using van't Hoff equation.

Note: Units are documented in descriptions but not enforced via DynamicQuantities
due to non-SI units (atm) used in atmospheric chemistry conventions.
All concentrations are in mol/L (M), pressures in atm, temperatures in K.
"""

# =============================================================================
# Equilibrium Constants at 298 K (Table 7.4)
# =============================================================================

# Water autoionization
const K_W_298 = 1.0e-14  # M^2
const DH_K_W = 13.35 * 4184  # J mol^-1 (converted from kcal mol^-1)

# CO2 system
const K_C1_298 = 4.3e-7  # M (CO2.H2O <-> H+ + HCO3-)
const K_C2_298 = 4.7e-11  # M (HCO3- <-> H+ + CO3^2-)
const DH_K_C1 = 1.83 * 4184  # J mol^-1
const DH_K_C2 = 3.55 * 4184  # J mol^-1

# SO2 system
const K_S1_298 = 1.3e-2  # M (SO2.H2O <-> H+ + HSO3-)
const K_S2_298 = 6.6e-8  # M (HSO3- <-> H+ + SO3^2-)
const DH_K_S1 = -4.16 * 4184  # J mol^-1 (negative!)
const DH_K_S2 = -2.23 * 4184  # J mol^-1 (negative!)

# NH3 system
const K_A1_298 = 1.7e-5  # M (NH3.H2O <-> NH4+ + OH-)
const DH_K_A1 = 8.65 * 4184  # J mol^-1

# HNO3 system
const K_N1_298 = 15.4  # M (HNO3 <-> H+ + NO3-)
# Strong acid, essentially complete dissociation

# H2O2 system
const K_H1_298 = 2.2e-12  # M (H2O2 <-> H+ + HO2-)
# Very weak, can usually be neglected

# =============================================================================
# Water Equilibrium
# =============================================================================

"""
    WaterEquilibrium(; name=:WaterEq)

Water autoionization equilibrium with temperature dependence.

Implements:
- Eq 7.10-7.12: H2O <-> H+ + OH-
- Eq 7.13: pH = -log10[H+]
- Table 7.4: K_w(T) = K_w(298) * exp((-dH/R) * (1/T - 1/298))

Variables:
- T: Temperature (K)
- H_plus: Hydrogen ion concentration (mol/L)
- OH_minus: Hydroxide ion concentration (mol/L)
- K_w: Water dissociation constant (mol^2/L^2)
- pH: Negative log of H+ concentration (dimensionless)
"""
@component function WaterEquilibrium(; name=:WaterEq)
    @constants begin
        R_gas = 8.314, [description = "Gas constant (J/mol/K)"]
        T_ref = 298.0, [description = "Reference temperature (K)"]
        K_w_298 = 1.0e-14, [description = "K_w at 298 K (M^2)"]
        dH_Kw = 55856.4, [description = "Enthalpy for K_w (13.35 kcal/mol = 55856 J/mol)"]
    end

    @variables begin
        T(t), [description = "Temperature (K)"]
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
        OH_minus(t), [description = "Hydroxide ion concentration (mol/L)"]
        K_w(t), [description = "Water dissociation constant (M^2)"]
        pH(t), [description = "pH value (dimensionless)"]
    end

    eqs = [
        # Temperature dependence of K_w (van't Hoff equation)
        K_w ~ K_w_298 * exp((-dH_Kw / R_gas) * (1/T - 1/T_ref)),

        # Eq 7.11: Water equilibrium
        K_w ~ H_plus * OH_minus,

        # Eq 7.13: pH definition
        pH ~ -log10(H_plus),
    ]

    return System(eqs, t; name)
end

# =============================================================================
# CO2 Equilibria (Section 7.3.2)
# =============================================================================

"""
    CO2Equilibria(; name=:CO2Eq)

CO2 dissolution and dissociation equilibria with temperature dependence.

Implements:
- Eq 7.14-7.16: CO2(g) <-> CO2.H2O <-> H+ + HCO3- <-> H+ + CO3^2-
- Eq 7.17: H_CO2 = [CO2.H2O] / p_CO2
- Eq 7.18: K_c1 = [H+][HCO3-] / [CO2.H2O]
- Eq 7.19: K_c2 = [H+][CO3^2-] / [HCO3-]
- Eq 7.20-7.22: Species concentrations
- Eq 7.24: H*_CO2 = H_CO2 * (1 + K_c1/[H+] + K_c1*K_c2/[H+]^2)

Variables:
- T: Temperature (K)
- H_plus: Hydrogen ion concentration (mol/L)
- p_CO2: Partial pressure of CO2 (atm)
- CO2_aq: [CO2.H2O] concentration (mol/L)
- HCO3_minus: Bicarbonate concentration (mol/L)
- CO3_2minus: Carbonate concentration (mol/L)
- C_total: Total dissolved inorganic carbon (mol/L)
- H_CO2_eff: Effective Henry's law constant (M/atm)
"""
@component function CO2Equilibria(; name=:CO2Eq)
    @constants begin
        R_gas = 8.314, [description = "Gas constant (J/mol/K)"]
        T_ref = 298.0, [description = "Reference temperature (K)"]
        H_CO2_298 = 3.4e-2, [description = "Henry's law constant for CO2 at 298 K (M/atm)"]
        dH_H_CO2 = -20292.4, [description = "Heat of dissolution for CO2 (-4.85 kcal/mol)"]
        K_c1_298 = 4.3e-7, [description = "First dissociation constant at 298 K (M)"]
        K_c2_298 = 4.7e-11, [description = "Second dissociation constant at 298 K (M)"]
        dH_Kc1 = 7657.12, [description = "Enthalpy for K_c1 (1.83 kcal/mol)"]
        dH_Kc2 = 14853.2, [description = "Enthalpy for K_c2 (3.55 kcal/mol)"]
    end

    @variables begin
        T(t), [description = "Temperature (K)"]
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
        p_CO2(t), [description = "Partial pressure of CO2 (atm)"]
        CO2_aq(t), [description = "Aqueous CO2 concentration [CO2.H2O] (mol/L)"]
        HCO3_minus(t), [description = "Bicarbonate concentration [HCO3-] (mol/L)"]
        CO3_2minus(t), [description = "Carbonate concentration [CO3^2-] (mol/L)"]
        C_total(t), [description = "Total dissolved inorganic carbon (mol/L)"]
        H_CO2(t), [description = "Henry's law constant for CO2 (M/atm)"]
        K_c1(t), [description = "First dissociation constant (M)"]
        K_c2(t), [description = "Second dissociation constant (M)"]
        H_CO2_eff(t), [description = "Effective Henry's law constant for CO2 (M/atm)"]
    end

    eqs = [
        # Temperature-dependent equilibrium constants
        H_CO2 ~ H_CO2_298 * exp((dH_H_CO2 / R_gas) * (1/T_ref - 1/T)),
        K_c1 ~ K_c1_298 * exp((-dH_Kc1 / R_gas) * (1/T - 1/T_ref)),
        K_c2 ~ K_c2_298 * exp((-dH_Kc2 / R_gas) * (1/T - 1/T_ref)),

        # Eq 7.17: Henry's law for CO2
        CO2_aq ~ H_CO2 * p_CO2,

        # Eq 7.18: First dissociation equilibrium
        # K_c1 = [H+][HCO3-] / [CO2.H2O]
        HCO3_minus ~ K_c1 * CO2_aq / H_plus,

        # Eq 7.19: Second dissociation equilibrium
        # K_c2 = [H+][CO3^2-] / [HCO3-]
        CO3_2minus ~ K_c2 * HCO3_minus / H_plus,

        # Eq 7.20-7.22: Total dissolved carbon
        C_total ~ CO2_aq + HCO3_minus + CO3_2minus,

        # Eq 7.24: Effective Henry's law constant
        H_CO2_eff ~ H_CO2 * (1 + K_c1/H_plus + K_c1*K_c2/H_plus^2),
    ]

    return System(eqs, t; name)
end

# =============================================================================
# SO2 Equilibria (Section 7.3.3)
# =============================================================================

"""
    SO2Equilibria(; name=:SO2Eq)

SO2 dissolution and dissociation equilibria with temperature dependence.
This is the S(IV) system.

Implements:
- Eq 7.30-7.32: SO2(g) <-> SO2.H2O <-> H+ + HSO3- <-> H+ + SO3^2-
- Eq 7.33: H_SO2 = [SO2.H2O] / p_SO2
- Eq 7.34: K_s1 = [H+][HSO3-] / [SO2.H2O]
- Eq 7.35: K_s2 = [H+][SO3^2-] / [HSO3-]
- Eq 7.36-7.38: S(IV) species concentrations
- Eq 7.39: [S(IV)] = [SO2.H2O] + [HSO3-] + [SO3^2-]
- Eq 7.40: [S(IV)] = H_SO2 * p_SO2 * (1 + K_s1/[H+] + K_s1*K_s2/[H+]^2)
- Eq 7.41: H*_S(IV) = H_SO2 * (1 + K_s1/[H+] + K_s1*K_s2/[H+]^2)
- Eq 7.43-7.45: S(IV) mole fractions (alpha values)

Variables:
- T: Temperature (K)
- H_plus: Hydrogen ion concentration (mol/L)
- p_SO2: Partial pressure of SO2 (atm)
- SO2_aq: [SO2.H2O] concentration (mol/L)
- HSO3_minus: Bisulfite concentration (mol/L)
- SO3_2minus: Sulfite concentration (mol/L)
- S_IV_total: Total S(IV) concentration (mol/L)
- H_SO2_eff: Effective Henry's law constant (M/atm)
- alpha_0: Mole fraction of SO2.H2O (dimensionless)
- alpha_1: Mole fraction of HSO3- (dimensionless)
- alpha_2: Mole fraction of SO3^2- (dimensionless)
"""
@component function SO2Equilibria(; name=:SO2Eq)
    @constants begin
        R_gas = 8.314, [description = "Gas constant (J/mol/K)"]
        T_ref = 298.0, [description = "Reference temperature (K)"]
        H_SO2_298 = 1.23, [description = "Henry's law constant for SO2 at 298 K (M/atm)"]
        dH_H_SO2 = -26150.0, [description = "Heat of dissolution for SO2 (-6.25 kcal/mol)"]
        K_s1_298 = 1.3e-2, [description = "First dissociation constant at 298 K (M)"]
        K_s2_298 = 6.6e-8, [description = "Second dissociation constant at 298 K (M)"]
        dH_Ks1 = -17405.44, [description = "Enthalpy for K_s1 (-4.16 kcal/mol)"]
        dH_Ks2 = -9330.32, [description = "Enthalpy for K_s2 (-2.23 kcal/mol)"]
    end

    @variables begin
        T(t), [description = "Temperature (K)"]
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
        p_SO2(t), [description = "Partial pressure of SO2 (atm)"]
        SO2_aq(t), [description = "Aqueous SO2 concentration [SO2.H2O] (mol/L)"]
        HSO3_minus(t), [description = "Bisulfite concentration [HSO3-] (mol/L)"]
        SO3_2minus(t), [description = "Sulfite concentration [SO3^2-] (mol/L)"]
        S_IV_total(t), [description = "Total S(IV) concentration (mol/L)"]
        H_SO2(t), [description = "Henry's law constant for SO2 (M/atm)"]
        K_s1(t), [description = "First dissociation constant (M)"]
        K_s2(t), [description = "Second dissociation constant (M)"]
        H_SO2_eff(t), [description = "Effective Henry's law constant for SO2 (M/atm)"]
        alpha_0(t), [description = "Mole fraction of SO2.H2O (dimensionless)"]
        alpha_1(t), [description = "Mole fraction of HSO3- (dimensionless)"]
        alpha_2(t), [description = "Mole fraction of SO3^2- (dimensionless)"]
    end

    eqs = [
        # Temperature-dependent equilibrium constants (note: negative dH means K increases with T)
        H_SO2 ~ H_SO2_298 * exp((dH_H_SO2 / R_gas) * (1/T_ref - 1/T)),
        K_s1 ~ K_s1_298 * exp((-dH_Ks1 / R_gas) * (1/T - 1/T_ref)),
        K_s2 ~ K_s2_298 * exp((-dH_Ks2 / R_gas) * (1/T - 1/T_ref)),

        # Eq 7.33: Henry's law for SO2
        SO2_aq ~ H_SO2 * p_SO2,

        # Eq 7.34: First dissociation equilibrium
        # K_s1 = [H+][HSO3-] / [SO2.H2O]
        HSO3_minus ~ K_s1 * SO2_aq / H_plus,

        # Eq 7.35: Second dissociation equilibrium
        # K_s2 = [H+][SO3^2-] / [HSO3-]
        SO3_2minus ~ K_s2 * HSO3_minus / H_plus,

        # Eq 7.39: Total S(IV)
        S_IV_total ~ SO2_aq + HSO3_minus + SO3_2minus,

        # Eq 7.41: Effective Henry's law constant
        H_SO2_eff ~ H_SO2 * (1 + K_s1/H_plus + K_s1*K_s2/H_plus^2),

        # Eq 7.43-7.45: Mole fractions (alpha values)
        # alpha_0 = [H+]^2 / ([H+]^2 + K_s1*[H+] + K_s1*K_s2)
        alpha_0 ~ H_plus^2 / (H_plus^2 + K_s1*H_plus + K_s1*K_s2),

        # alpha_1 = K_s1*[H+] / ([H+]^2 + K_s1*[H+] + K_s1*K_s2)
        alpha_1 ~ K_s1*H_plus / (H_plus^2 + K_s1*H_plus + K_s1*K_s2),

        # alpha_2 = K_s1*K_s2 / ([H+]^2 + K_s1*[H+] + K_s1*K_s2)
        alpha_2 ~ K_s1*K_s2 / (H_plus^2 + K_s1*H_plus + K_s1*K_s2),
    ]

    return System(eqs, t; name)
end

# =============================================================================
# NH3 Equilibria (Section 7.3.4)
# =============================================================================

"""
    NH3Equilibria(; name=:NH3Eq)

NH3 dissolution and dissociation equilibria with temperature dependence.

Implements:
- Eq 7.46-7.47: NH3(g) <-> NH3.H2O <-> NH4+ + OH-
- Eq 7.48: H_NH3 = [NH3.H2O] / p_NH3
- Eq 7.49: K_a1 = [NH4+][OH-] / [NH3.H2O]
- Eq 7.50: [NH4+] = (H_NH3 * K_a1 / K_w) * p_NH3 * [H+]
- Eq 7.51: [NH3_T] = H_NH3 * p_NH3 * (1 + K_a1*[H+]/K_w)

Variables:
- T: Temperature (K)
- H_plus: Hydrogen ion concentration (mol/L)
- OH_minus: Hydroxide ion concentration (mol/L)
- K_w: Water dissociation constant (M^2)
- p_NH3: Partial pressure of NH3 (atm)
- NH3_aq: [NH3.H2O] concentration (mol/L)
- NH4_plus: Ammonium concentration (mol/L)
- NH3_total: Total dissolved ammonia (mol/L)
- H_NH3_eff: Effective Henry's law constant (M/atm)
"""
@component function NH3Equilibria(; name=:NH3Eq)
    @constants begin
        R_gas = 8.314, [description = "Gas constant (J/mol/K)"]
        T_ref = 298.0, [description = "Reference temperature (K)"]
        H_NH3_298 = 62.0, [description = "Henry's law constant for NH3 at 298 K (M/atm)"]
        dH_H_NH3 = -34183.28, [description = "Heat of dissolution for NH3 (-8.17 kcal/mol)"]
        K_a1_298 = 1.7e-5, [description = "Dissociation constant at 298 K (M)"]
        dH_Ka1 = 36191.6, [description = "Enthalpy for K_a1 (8.65 kcal/mol)"]
        K_w_298 = 1.0e-14, [description = "K_w at 298 K (M^2)"]
        dH_Kw = 55856.4, [description = "Enthalpy for K_w (13.35 kcal/mol)"]
    end

    @variables begin
        T(t), [description = "Temperature (K)"]
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
        OH_minus(t), [description = "Hydroxide ion concentration (mol/L)"]
        K_w(t), [description = "Water dissociation constant (M^2)"]
        p_NH3(t), [description = "Partial pressure of NH3 (atm)"]
        NH3_aq(t), [description = "Aqueous NH3 concentration [NH3.H2O] (mol/L)"]
        NH4_plus(t), [description = "Ammonium concentration [NH4+] (mol/L)"]
        NH3_total(t), [description = "Total dissolved ammonia (mol/L)"]
        H_NH3(t), [description = "Henry's law constant for NH3 (M/atm)"]
        K_a1(t), [description = "Dissociation constant (M)"]
        H_NH3_eff(t), [description = "Effective Henry's law constant for NH3 (M/atm)"]
    end

    eqs = [
        # Temperature-dependent equilibrium constants
        H_NH3 ~ H_NH3_298 * exp((dH_H_NH3 / R_gas) * (1/T_ref - 1/T)),
        K_a1 ~ K_a1_298 * exp((-dH_Ka1 / R_gas) * (1/T - 1/T_ref)),
        K_w ~ K_w_298 * exp((-dH_Kw / R_gas) * (1/T - 1/T_ref)),

        # Water equilibrium
        OH_minus ~ K_w / H_plus,

        # Eq 7.48: Henry's law for NH3
        NH3_aq ~ H_NH3 * p_NH3,

        # Eq 7.49: Dissociation equilibrium
        # K_a1 = [NH4+][OH-] / [NH3.H2O]
        # Rearranged: [NH4+] = K_a1 * [NH3.H2O] / [OH-]
        NH4_plus ~ K_a1 * NH3_aq / OH_minus,

        # Eq 7.51: Total dissolved ammonia
        # [NH3_T] = [NH3.H2O] + [NH4+] = H_NH3 * p_NH3 * (1 + K_a1*[H+]/K_w)
        NH3_total ~ NH3_aq + NH4_plus,

        # Effective Henry's law constant
        # H*_NH3 = H_NH3 * (1 + K_a1*[H+]/K_w)
        H_NH3_eff ~ H_NH3 * (1 + K_a1*H_plus/K_w),
    ]

    return System(eqs, t; name)
end

# =============================================================================
# HNO3 Equilibria (Section 7.3.5)
# =============================================================================

"""
    HNO3Equilibria(; name=:HNO3Eq)

HNO3 dissolution and dissociation equilibria.

Implements:
- Eq 7.53-7.54: HNO3(g) <-> HNO3(aq) <-> H+ + NO3-
- H_HNO3 = [HNO3(aq)] / p_HNO3 = 2.1e5 M/atm at 298 K
- K_n1 = [H+][NO3-] / [HNO3(aq)] = 15.4 M at 298 K
- Eq 7.58: [NO3-] ~ H_HNO3 * K_n1 / [H+] * p_HNO3 (for complete dissociation)
- Eq 7.61: H*_HNO3 ~ H_HNO3 * K_n1 / [H+] = 3.2e6 / [H+]

Note: HNO3 is a strong acid and essentially completely dissociates in water.

Variables:
- T: Temperature (K)
- H_plus: Hydrogen ion concentration (mol/L)
- p_HNO3: Partial pressure of HNO3 (atm)
- HNO3_aq: Undissociated HNO3 concentration (mol/L)
- NO3_minus: Nitrate concentration (mol/L)
- HNO3_total: Total dissolved nitric acid (mol/L)
- H_HNO3_eff: Effective Henry's law constant (M/atm)
"""
@component function HNO3Equilibria(; name=:HNO3Eq)
    @constants begin
        R_gas = 8.314, [description = "Gas constant (J/mol/K)"]
        T_ref = 298.0, [description = "Reference temperature (K)"]
        H_HNO3_298 = 2.1e5, [description = "Henry's law constant for HNO3 at 298 K (M/atm)"]
        K_n1_298 = 15.4, [description = "Dissociation constant at 298 K (M)"]
    end

    @variables begin
        T(t), [description = "Temperature (K)"]
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
        p_HNO3(t), [description = "Partial pressure of HNO3 (atm)"]
        HNO3_aq(t), [description = "Undissociated HNO3 concentration (mol/L)"]
        NO3_minus(t), [description = "Nitrate concentration [NO3-] (mol/L)"]
        HNO3_total(t), [description = "Total dissolved nitric acid (mol/L)"]
        H_HNO3(t), [description = "Henry's law constant for HNO3 (M/atm)"]
        K_n1(t), [description = "Dissociation constant (M)"]
        H_HNO3_eff(t), [description = "Effective Henry's law constant for HNO3 (M/atm)"]
    end

    eqs = [
        # Temperature-dependent constants (assuming weak temperature dependence for HNO3)
        H_HNO3 ~ H_HNO3_298,
        K_n1 ~ K_n1_298,

        # Henry's law for undissociated HNO3
        HNO3_aq ~ H_HNO3 * p_HNO3,

        # Dissociation equilibrium
        # K_n1 = [H+][NO3-] / [HNO3(aq)]
        NO3_minus ~ K_n1 * HNO3_aq / H_plus,

        # Total dissolved nitric acid
        HNO3_total ~ HNO3_aq + NO3_minus,

        # Eq 7.61: Effective Henry's law constant
        # H*_HNO3 ~ H_HNO3 * (1 + K_n1/[H+])
        H_HNO3_eff ~ H_HNO3 * (1 + K_n1/H_plus),
    ]

    return System(eqs, t; name)
end

# =============================================================================
# H2O2 Equilibria (Section 7.3.6)
# =============================================================================

"""
    H2O2Equilibria(; name=:H2O2Eq)

H2O2 dissolution and (weak) dissociation equilibria.

Implements:
- H_H2O2 = [H2O2(aq)] / p_H2O2 = 1e5 M/atm at 298 K
- K_h1 = [H+][HO2-] / [H2O2(aq)] = 2.2e-12 M (very weak, usually neglected)

Note: H2O2 is a very weak acid and dissociation is typically neglected.

Variables:
- T: Temperature (K)
- H_plus: Hydrogen ion concentration (mol/L)
- p_H2O2: Partial pressure of H2O2 (atm)
- H2O2_aq: Aqueous H2O2 concentration (mol/L)
- HO2_minus: Hydroperoxide anion concentration (mol/L)
- H2O2_total: Total dissolved hydrogen peroxide (mol/L)
- H_H2O2_eff: Effective Henry's law constant (M/atm)
"""
@component function H2O2Equilibria(; name=:H2O2Eq)
    @constants begin
        R_gas = 8.314, [description = "Gas constant (J/mol/K)"]
        T_ref = 298.0, [description = "Reference temperature (K)"]
        H_H2O2_298 = 1.0e5, [description = "Henry's law constant for H2O2 at 298 K (M/atm)"]
        dH_H_H2O2 = -60668.0, [description = "Heat of dissolution for H2O2 (-14.5 kcal/mol)"]
        K_h1_298 = 2.2e-12, [description = "Dissociation constant at 298 K (very weak) (M)"]
    end

    @variables begin
        T(t), [description = "Temperature (K)"]
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
        p_H2O2(t), [description = "Partial pressure of H2O2 (atm)"]
        H2O2_aq(t), [description = "Aqueous H2O2 concentration (mol/L)"]
        HO2_minus(t), [description = "Hydroperoxide anion concentration [HO2-] (mol/L)"]
        H2O2_total(t), [description = "Total dissolved hydrogen peroxide (mol/L)"]
        H_H2O2(t), [description = "Henry's law constant for H2O2 (M/atm)"]
        K_h1(t), [description = "Dissociation constant (very weak) (M)"]
        H_H2O2_eff(t), [description = "Effective Henry's law constant for H2O2 (M/atm)"]
    end

    eqs = [
        # Temperature-dependent Henry's law constant
        H_H2O2 ~ H_H2O2_298 * exp((dH_H_H2O2 / R_gas) * (1/T_ref - 1/T)),
        K_h1 ~ K_h1_298,

        # Henry's law for H2O2
        H2O2_aq ~ H_H2O2 * p_H2O2,

        # Dissociation equilibrium (very weak, often negligible)
        # K_h1 = [H+][HO2-] / [H2O2(aq)]
        HO2_minus ~ K_h1 * H2O2_aq / H_plus,

        # Total dissolved H2O2
        H2O2_total ~ H2O2_aq + HO2_minus,

        # Effective Henry's law constant
        H_H2O2_eff ~ H_H2O2 * (1 + K_h1/H_plus),
    ]

    return System(eqs, t; name)
end

# =============================================================================
# O3 Equilibria (for oxidation reactions)
# =============================================================================

"""
    O3Equilibria(; name=:O3Eq)

O3 dissolution equilibrium with temperature dependence.

O3 does not dissociate in water, so only Henry's law applies.

Variables:
- T: Temperature (K)
- p_O3: Partial pressure of O3 (atm)
- O3_aq: Aqueous O3 concentration (mol/L)
- H_O3: Henry's law constant (M/atm)
"""
@component function O3Equilibria(; name=:O3Eq)
    @constants begin
        R_gas = 8.314, [description = "Gas constant (J/mol/K)"]
        T_ref = 298.0, [description = "Reference temperature (K)"]
        H_O3_298 = 1.1e-2, [description = "Henry's law constant for O3 at 298 K (M/atm)"]
        dH_H_O3 = -21087.36, [description = "Heat of dissolution for O3 (-5.04 kcal/mol)"]
    end

    @variables begin
        T(t), [description = "Temperature (K)"]
        p_O3(t), [description = "Partial pressure of O3 (atm)"]
        O3_aq(t), [description = "Aqueous O3 concentration (mol/L)"]
        H_O3(t), [description = "Henry's law constant for O3 (M/atm)"]
    end

    eqs = [
        # Temperature-dependent Henry's law constant
        H_O3 ~ H_O3_298 * exp((dH_H_O3 / R_gas) * (1/T_ref - 1/T)),

        # Henry's law for O3
        O3_aq ~ H_O3 * p_O3,
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Combined Aqueous Equilibria System
# =============================================================================

"""
    AqueousEquilibria(; name=:AqEquilibria)

Combined aqueous equilibria system for atmospheric cloud droplet chemistry.

This component combines all individual equilibrium systems into a single
comprehensive system that can be used for cloud chemistry calculations.

Subsystems:
- water_eq: Water autoionization
- co2_eq: CO2/bicarbonate/carbonate equilibria
- so2_eq: SO2/bisulfite/sulfite (S(IV)) equilibria
- nh3_eq: NH3/ammonium equilibria
- hno3_eq: HNO3/nitrate equilibria
- h2o2_eq: H2O2 equilibria
- o3_eq: O3 dissolution
"""
@component function AqueousEquilibria(; name=:AqEquilibria)
    # Create subsystems
    water_eq = WaterEquilibrium(; name=:water_eq)
    co2_eq = CO2Equilibria(; name=:co2_eq)
    so2_eq = SO2Equilibria(; name=:so2_eq)
    nh3_eq = NH3Equilibria(; name=:nh3_eq)
    hno3_eq = HNO3Equilibria(; name=:hno3_eq)
    h2o2_eq = H2O2Equilibria(; name=:h2o2_eq)
    o3_eq = O3Equilibria(; name=:o3_eq)

    @variables begin
        T(t), [description = "Temperature (K)"]
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
    end

    # Connection equations: share T and H_plus across all subsystems
    eqs = [
        # Temperature coupling
        water_eq.T ~ T,
        co2_eq.T ~ T,
        so2_eq.T ~ T,
        nh3_eq.T ~ T,
        hno3_eq.T ~ T,
        h2o2_eq.T ~ T,
        o3_eq.T ~ T,

        # H+ coupling (all subsystems use the same H+ from water equilibrium)
        water_eq.H_plus ~ H_plus,
        co2_eq.H_plus ~ H_plus,
        so2_eq.H_plus ~ H_plus,
        nh3_eq.H_plus ~ H_plus,
        hno3_eq.H_plus ~ H_plus,
        h2o2_eq.H_plus ~ H_plus,
    ]

    return System(eqs, t;
        systems=[water_eq, co2_eq, so2_eq, nh3_eq, hno3_eq, h2o2_eq, o3_eq],
        name)
end
