"""
    S(IV) Oxidation Kinetics (Sulfate Formation)

Implements S(IV) to S(VI) oxidation reaction kinetics from Seinfeld & Pandis
Chapter 7, Section 7.5.

Key oxidation pathways:
1. S(IV) + O3 (Eq 7.79-7.80)
2. S(IV) + H2O2 (Eq 7.84)
3. S(IV) + O2 catalyzed by Fe(III) (Eq 7.92-7.96)
4. S(IV) + O2 catalyzed by Mn(II) (Eq 7.98-7.101)
5. Fe/Mn synergism (Eq 7.102)

Also includes rate conversion equations (Section 7.4):
- Eq 7.75: R''_a (ppb h^-1) = 3.6e6 * L * R * T * R_a (M s^-1)
- Eq 7.76: R_bar_a (%h^-1) = 3.6e8 * R_a * L * R * T / xi_SO2
- Eq 7.78: tau_SO2 (s) = xi_SO2 / (1e3 * R_a * L * R * T)

Note: Units are documented in descriptions but not enforced via DynamicQuantities
due to non-SI units (atm) used in atmospheric chemistry conventions.
"""

# =============================================================================
# Rate Constants for S(IV) Oxidation
# =============================================================================

# O3 oxidation rate constants (Eq 7.79)
const K0_O3 = 2.4e4   # M^-1 s^-1 for SO2.H2O + O3
const K1_O3 = 3.7e5   # M^-1 s^-1 for HSO3- + O3
const K2_O3 = 1.5e9   # M^-1 s^-1 for SO3^2- + O3

# H2O2 oxidation rate constants (Eq 7.84)
const K_H2O2 = 7.5e7  # M^-2 s^-1
const K_H2O2_eq = 13.0  # M^-1 (equilibrium constant in denominator)

# Fe(III) catalyzed oxidation rate constants
const K_FE_LOW_PH = 6.0  # s^-1 for pH 0-3.6
const K_FE_PH4 = 1.0e9   # M^-1 s^-1 for pH ~4.0
const K_FE_PH5_6 = 1.0e-3  # s^-1 for pH 5.0-6.0
const K_FE_PH7 = 1.0e-4  # s^-1 for pH 7.0

# Mn(II) catalyzed oxidation rate constants
const K0_MN = 1000.0  # M^-1 s^-1

# Fe/Mn synergism rate constants (Eq 7.102)
const K_MN_SIV = 750.0      # M^-1 s^-1
const K_FE_SIV = 2600.0     # M^-1 s^-1
const K_FEMN_SIV = 1.0e10   # M^-2 s^-1

# =============================================================================
# S(IV) + O3 Oxidation (Section 7.5.1)
# =============================================================================

"""
    SulfateFormationO3(; name=:SulfateO3)

S(IV) oxidation by O3 in aqueous phase.

Implements Eq 7.79-7.80:
R_O3 = (k0*[SO2.H2O] + k1*[HSO3-] + k2*[SO3^2-]) * [O3(aq)]

where:
- k0 = 2.4e4 M^-1 s^-1 (SO2.H2O + O3)
- k1 = 3.7e5 M^-1 s^-1 (HSO3- + O3)
- k2 = 1.5e9 M^-1 s^-1 (SO3^2- + O3)

The reaction is fast at high pH due to the large k2 value for SO3^2-.

Variables:
- SO2_aq: [SO2.H2O] concentration (mol/L)
- HSO3_minus: [HSO3-] concentration (mol/L)
- SO3_2minus: [SO3^2-] concentration (mol/L)
- O3_aq: [O3(aq)] concentration (mol/L)
- R_O3: Oxidation rate (mol/L/s)
"""
@component function SulfateFormationO3(; name=:SulfateO3)
    @constants begin
        k0 = 2.4e4, [description = "Rate constant for SO2.H2O + O3 (M^-1 s^-1)"]
        k1 = 3.7e5, [description = "Rate constant for HSO3- + O3 (M^-1 s^-1)"]
        k2 = 1.5e9, [description = "Rate constant for SO3^2- + O3 (M^-1 s^-1)"]
    end

    @variables begin
        SO2_aq(t), [description = "Aqueous SO2 concentration [SO2.H2O] (mol/L)"]
        HSO3_minus(t), [description = "Bisulfite concentration [HSO3-] (mol/L)"]
        SO3_2minus(t), [description = "Sulfite concentration [SO3^2-] (mol/L)"]
        O3_aq(t), [description = "Aqueous O3 concentration (mol/L)"]
        R_O3(t), [description = "S(IV) oxidation rate by O3 (mol/L/s)"]
    end

    eqs = [
        # Eq 7.79: Total oxidation rate by O3
        # R = (k0*[SO2.H2O] + k1*[HSO3-] + k2*[SO3^2-]) * [O3(aq)]
        R_O3 ~ (k0*SO2_aq + k1*HSO3_minus + k2*SO3_2minus) * O3_aq,
    ]

    return System(eqs, t; name)
end

# =============================================================================
# S(IV) + H2O2 Oxidation (Section 7.5.2)
# =============================================================================

"""
    SulfateFormationH2O2(; name=:SulfateH2O2)

S(IV) oxidation by H2O2 in aqueous phase.

Implements Eq 7.84:
R_H2O2 = k*[H+]*[H2O2]*[HSO3-] / (1 + K*[H+])

where:
- k = 7.5e7 M^-2 s^-1
- K = 13 M^-1

This reaction is relatively pH-independent due to the [H+] terms canceling,
making it important across a wide pH range.

Variables:
- H_plus: Hydrogen ion concentration (mol/L)
- H2O2_aq: Aqueous H2O2 concentration (mol/L)
- HSO3_minus: Bisulfite concentration (mol/L)
- R_H2O2: Oxidation rate (mol/L/s)
"""
@component function SulfateFormationH2O2(; name=:SulfateH2O2)
    @constants begin
        k = 7.5e7, [description = "Rate constant for H2O2 oxidation (M^-2 s^-1)"]
        K_eq = 13.0, [description = "Equilibrium constant in denominator (M^-1)"]
    end

    @variables begin
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
        H2O2_aq(t), [description = "Aqueous H2O2 concentration (mol/L)"]
        HSO3_minus(t), [description = "Bisulfite concentration [HSO3-] (mol/L)"]
        R_H2O2(t), [description = "S(IV) oxidation rate by H2O2 (mol/L/s)"]
    end

    eqs = [
        # Eq 7.84: H2O2 oxidation rate
        # R = k*[H+]*[H2O2]*[HSO3-] / (1 + K*[H+])
        R_H2O2 ~ k * H_plus * H2O2_aq * HSO3_minus / (1 + K_eq * H_plus),
    ]

    return System(eqs, t; name)
end

# =============================================================================
# S(IV) + O2 Catalyzed by Fe(III) (Section 7.5.3)
# =============================================================================

"""
    SulfateFormationFe(; name=:SulfateFe)

S(IV) oxidation by O2 catalyzed by Fe(III) in aqueous phase.

The rate expression depends strongly on pH:

For pH 0-3.6 (Eq 7.92):
R_Fe = k*[Fe(III)]*[S(IV)] / [H+], k = 6 s^-1

For pH ~4.0 (Eq 7.93):
R_Fe = 1e9 * [S(IV)] * [Fe(III)]^2

For pH 5.0-6.0 (Eq 7.94):
R_Fe = 1e-3 * [S(IV)]

For pH 7.0 (Eq 7.95):
R_Fe = 1e-4 * [S(IV)]

This implementation uses a simplified approach with the low-pH rate law,
which is most relevant for typical cloud droplet conditions.

Variables:
- H_plus: Hydrogen ion concentration (mol/L)
- Fe_III: Fe(III) concentration (mol/L)
- S_IV_total: Total S(IV) concentration (mol/L)
- R_Fe: Oxidation rate (mol/L/s)
"""
@component function SulfateFormationFe(; name=:SulfateFe)
    @constants begin
        k_low_pH = 6.0, [description = "Rate constant for pH 0-3.6 (s^-1)"]
    end

    @variables begin
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
        Fe_III(t), [description = "Fe(III) concentration (mol/L)"]
        S_IV_total(t), [description = "Total S(IV) concentration (mol/L)"]
        R_Fe(t), [description = "S(IV) oxidation rate by Fe(III)-catalyzed O2 (mol/L/s)"]
    end

    eqs = [
        # Eq 7.92: Fe(III)-catalyzed oxidation (low pH regime)
        # R = k*[Fe(III)]*[S(IV)] / [H+]
        R_Fe ~ k_low_pH * Fe_III * S_IV_total / H_plus,
    ]

    return System(eqs, t; name)
end

# =============================================================================
# S(IV) + O2 Catalyzed by Mn(II) (Section 7.5.4)
# =============================================================================

"""
    SulfateFormationMn(; name=:SulfateMn)

S(IV) oxidation by O2 catalyzed by Mn(II) in aqueous phase.

The rate expression depends on S(IV) concentration:

For high S(IV) (>100 uM) (Eq 7.98):
R_Mn = k0 * [Mn]^2  (zero-order in S(IV))

For low S(IV) (<1 uM) (Eq 7.101):
R_Mn = k0 * [Mn] * [S(IV)], k0 = 1000 M^-1 s^-1

This implementation uses the low S(IV) rate law, which is first-order in
both Mn and S(IV).

Variables:
- Mn_II: Mn(II) concentration (mol/L)
- S_IV_total: Total S(IV) concentration (mol/L)
- R_Mn: Oxidation rate (mol/L/s)
"""
@component function SulfateFormationMn(; name=:SulfateMn)
    @constants begin
        k0 = 1000.0, [description = "Rate constant for low S(IV) (M^-1 s^-1)"]
    end

    @variables begin
        Mn_II(t), [description = "Mn(II) concentration (mol/L)"]
        S_IV_total(t), [description = "Total S(IV) concentration (mol/L)"]
        R_Mn(t), [description = "S(IV) oxidation rate by Mn(II)-catalyzed O2 (mol/L/s)"]
    end

    eqs = [
        # Eq 7.101: Mn(II)-catalyzed oxidation (low S(IV) regime)
        # R = k0 * [Mn] * [S(IV)]
        R_Mn ~ k0 * Mn_II * S_IV_total,
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Fe/Mn Synergism (Section 7.5.5)
# =============================================================================

"""
    SulfateFormationFeMn(; name=:SulfateFeMn)

S(IV) oxidation by O2 with Fe/Mn synergism.

Implements Eq 7.102:
R = 750*[Mn(II)]*[S(IV)] + 2600*[Fe(III)]*[S(IV)] + 1.0e10*[Mn(II)]*[Fe(III)]*[S(IV)]

The synergistic term (10^10 coefficient) represents the greatly enhanced
oxidation rate when both Fe(III) and Mn(II) are present.

Variables:
- Mn_II: Mn(II) concentration (mol/L)
- Fe_III: Fe(III) concentration (mol/L)
- S_IV_total: Total S(IV) concentration (mol/L)
- R_Mn_term: Mn-only contribution (mol/L/s)
- R_Fe_term: Fe-only contribution (mol/L/s)
- R_synergy_term: Synergistic contribution (mol/L/s)
- R_FeMn: Total oxidation rate (mol/L/s)
"""
@component function SulfateFormationFeMn(; name=:SulfateFeMn)
    @constants begin
        k_Mn = 750.0, [description = "Rate constant for Mn-only term (M^-1 s^-1)"]
        k_Fe = 2600.0, [description = "Rate constant for Fe-only term (M^-1 s^-1)"]
        k_FeMn = 1.0e10, [description = "Rate constant for synergistic term (M^-2 s^-1)"]
    end

    @variables begin
        Mn_II(t), [description = "Mn(II) concentration (mol/L)"]
        Fe_III(t), [description = "Fe(III) concentration (mol/L)"]
        S_IV_total(t), [description = "Total S(IV) concentration (mol/L)"]
        R_Mn_term(t), [description = "Mn-only contribution to rate (mol/L/s)"]
        R_Fe_term(t), [description = "Fe-only contribution to rate (mol/L/s)"]
        R_synergy_term(t), [description = "Synergistic contribution to rate (mol/L/s)"]
        R_FeMn(t), [description = "Total S(IV) oxidation rate with Fe/Mn synergism (mol/L/s)"]
    end

    eqs = [
        # Eq 7.102: Three-term rate expression
        R_Mn_term ~ k_Mn * Mn_II * S_IV_total,
        R_Fe_term ~ k_Fe * Fe_III * S_IV_total,
        R_synergy_term ~ k_FeMn * Mn_II * Fe_III * S_IV_total,

        # Total rate
        R_FeMn ~ R_Mn_term + R_Fe_term + R_synergy_term,
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Combined Sulfate Formation System
# =============================================================================

"""
    SulfateFormation(; name=:SulfateForm)

Combined S(IV) to S(VI) oxidation kinetics including all pathways.

Subsystems:
- o3_ox: O3 oxidation pathway
- h2o2_ox: H2O2 oxidation pathway
- fe_ox: Fe(III)-catalyzed oxidation
- mn_ox: Mn(II)-catalyzed oxidation
- femn_ox: Fe/Mn synergistic oxidation

Variables:
- R_total: Total S(IV) oxidation rate (mol/L/s)
- S_VI_production: Sulfate production rate (mol/L/s)

Also calculates rate conversion factors (Section 7.4):
- R_ppb_hr: Rate in ppb hr^-1 (Eq 7.75)
- R_percent_hr: Rate in % hr^-1 (Eq 7.76)
- tau_SO2: SO2 lifetime (s) (Eq 7.78)
"""
@component function SulfateFormation(; name=:SulfateForm)
    # Create subsystems
    o3_ox = SulfateFormationO3(; name=:o3_ox)
    h2o2_ox = SulfateFormationH2O2(; name=:h2o2_ox)
    fe_ox = SulfateFormationFe(; name=:fe_ox)
    mn_ox = SulfateFormationMn(; name=:mn_ox)
    femn_ox = SulfateFormationFeMn(; name=:femn_ox)

    @constants begin
        R_atm = 0.08205, [description = "Gas constant for rate conversion (L*atm/mol/K)"]
    end

    @parameters begin
        L, [description = "Liquid water content (g/m^3)"]
        xi_SO2, [description = "SO2 mixing ratio (dimensionless)"]
    end

    @variables begin
        T(t), [description = "Temperature (K)"]
        H_plus(t), [description = "Hydrogen ion concentration (mol/L)"]
        SO2_aq(t), [description = "Aqueous SO2 concentration [SO2.H2O] (mol/L)"]
        HSO3_minus(t), [description = "Bisulfite concentration [HSO3-] (mol/L)"]
        SO3_2minus(t), [description = "Sulfite concentration [SO3^2-] (mol/L)"]
        S_IV_total(t), [description = "Total S(IV) concentration (mol/L)"]
        O3_aq(t), [description = "Aqueous O3 concentration (mol/L)"]
        H2O2_aq(t), [description = "Aqueous H2O2 concentration (mol/L)"]
        Fe_III(t), [description = "Fe(III) concentration (mol/L)"]
        Mn_II(t), [description = "Mn(II) concentration (mol/L)"]
        R_total(t), [description = "Total S(IV) oxidation rate (mol/L/s)"]
        S_VI_production(t), [description = "Sulfate production rate (mol/L/s)"]
        R_ppb_hr(t), [description = "Rate in ppb/hr (Eq 7.75)"]
        R_percent_hr(t), [description = "Rate in %/hr (Eq 7.76)"]
        tau_SO2(t), [description = "SO2 lifetime (s) (Eq 7.78)"]
    end

    # Connection equations
    eqs = [
        # Connect S(IV) species to O3 oxidation subsystem
        o3_ox.SO2_aq ~ SO2_aq,
        o3_ox.HSO3_minus ~ HSO3_minus,
        o3_ox.SO3_2minus ~ SO3_2minus,
        o3_ox.O3_aq ~ O3_aq,

        # Connect to H2O2 oxidation subsystem
        h2o2_ox.H_plus ~ H_plus,
        h2o2_ox.H2O2_aq ~ H2O2_aq,
        h2o2_ox.HSO3_minus ~ HSO3_minus,

        # Connect to Fe oxidation subsystem
        fe_ox.H_plus ~ H_plus,
        fe_ox.Fe_III ~ Fe_III,
        fe_ox.S_IV_total ~ S_IV_total,

        # Connect to Mn oxidation subsystem
        mn_ox.Mn_II ~ Mn_II,
        mn_ox.S_IV_total ~ S_IV_total,

        # Connect to Fe/Mn synergism subsystem
        femn_ox.Mn_II ~ Mn_II,
        femn_ox.Fe_III ~ Fe_III,
        femn_ox.S_IV_total ~ S_IV_total,

        # Total oxidation rate (sum of all pathways)
        # Note: femn_ox includes Fe and Mn terms, so we use it instead of fe_ox and mn_ox
        # to avoid double-counting when synergism is important
        R_total ~ o3_ox.R_O3 + h2o2_ox.R_H2O2 + femn_ox.R_FeMn,

        # Sulfate production equals S(IV) oxidation rate
        S_VI_production ~ R_total,

        # Eq 7.75: Rate conversion to ppb/hr
        # R''_a (ppb h^-1) = 3.6e6 * L * R * T * R_a (M s^-1)
        R_ppb_hr ~ 3.6e6 * L * R_atm * T * R_total,

        # Eq 7.76: Rate conversion to %/hr
        # R_bar_a (%h^-1) = 3.6e8 * R_a * L * R * T / xi_SO2
        R_percent_hr ~ 3.6e8 * R_total * L * R_atm * T / xi_SO2,

        # Eq 7.78: SO2 lifetime
        # tau_SO2 (s) = xi_SO2 / (1e3 * R_a * L * R * T)
        tau_SO2 ~ xi_SO2 / (1e3 * R_total * L * R_atm * T),
    ]

    return System(eqs, t;
        systems=[o3_ox, h2o2_ox, fe_ox, mn_ox, femn_ox],
        name)
end

# =============================================================================
# Rate Conversion Utilities (Section 7.4)
# =============================================================================

"""
    rate_to_ppb_hr(R_a, L, T)

Convert aqueous-phase rate to atmospheric mixing ratio rate (Eq 7.75).

Arguments:
- R_a: Aqueous reaction rate (M s^-1)
- L: Liquid water content (g m^-3)
- T: Temperature (K)

Returns:
- Rate in ppb hr^-1
"""
function rate_to_ppb_hr(R_a, L, T)
    R = 0.08205  # atm L mol^-1 K^-1
    return 3.6e6 * L * R * T * R_a
end

"""
    rate_to_percent_hr(R_a, L, T, xi_SO2)

Convert aqueous-phase rate to percent conversion rate (Eq 7.76).

Arguments:
- R_a: Aqueous reaction rate (M s^-1)
- L: Liquid water content (g m^-3)
- T: Temperature (K)
- xi_SO2: SO2 mixing ratio (dimensionless)

Returns:
- Rate in % hr^-1
"""
function rate_to_percent_hr(R_a, L, T, xi_SO2)
    R = 0.08205  # atm L mol^-1 K^-1
    return 3.6e8 * R_a * L * R * T / xi_SO2
end

"""
    so2_lifetime(R_a, L, T, xi_SO2)

Calculate SO2 lifetime in seconds (Eq 7.78).

Arguments:
- R_a: Aqueous reaction rate (M s^-1)
- L: Liquid water content (g m^-3)
- T: Temperature (K)
- xi_SO2: SO2 mixing ratio (dimensionless)

Returns:
- SO2 lifetime in seconds
"""
function so2_lifetime(R_a, L, T, xi_SO2)
    R = 0.08205  # atm L mol^-1 K^-1
    return xi_SO2 / (1e3 * R_a * L * R * T)
end
