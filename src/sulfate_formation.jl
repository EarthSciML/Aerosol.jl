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

All concentrations are in mol/m^3, pressures in Pa (SI), temperatures in K.
Rate conversion equations use R in J/(mol·K) (= Pa·m^3/(mol·K)).
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

  - SO2_aq: [SO2.H2O] concentration (mol/m^3)
  - HSO3_minus: [HSO3-] concentration (mol/m^3)
  - SO3_2minus: [SO3^2-] concentration (mol/m^3)
  - O3_aq: [O3(aq)] concentration (mol/m^3)
  - R_O3: Oxidation rate (mol/m^3/s)
"""
@component function SulfateFormationO3(; name = :SulfateO3)
    @constants begin
        k0 = 2.4e4 / 1000.0,
        [description = "Rate constant for SO2.H2O + O3", unit = u"m^3/mol/s"]
        k1 = 3.7e5 / 1000.0,
        [description = "Rate constant for HSO3- + O3", unit = u"m^3/mol/s"]
        k2 = 1.5e9 / 1000.0,
        [description = "Rate constant for SO3^2- + O3", unit = u"m^3/mol/s"]
    end

    @variables begin
        SO2_aq(t), [description = "Aqueous SO2 concentration [SO2.H2O]", unit = u"mol/m^3"]
        HSO3_minus(t), [description = "Bisulfite concentration [HSO3-]", unit = u"mol/m^3"]
        SO3_2minus(t), [description = "Sulfite concentration [SO3^2-]", unit = u"mol/m^3"]
        O3_aq(t), [description = "Aqueous O3 concentration", unit = u"mol/m^3"]
        R_O3(t), [description = "S(IV) oxidation rate by O3", unit = u"mol/m^3/s"]
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

  - H_plus: Hydrogen ion concentration (mol/m^3)
  - H2O2_aq: Aqueous H2O2 concentration (mol/m^3)
  - HSO3_minus: Bisulfite concentration (mol/m^3)
  - R_H2O2: Oxidation rate (mol/m^3/s)
"""
@component function SulfateFormationH2O2(; name = :SulfateH2O2)
    @constants begin
        k = 7.5e7 / 1e6,
        [description = "Rate constant for H2O2 oxidation", unit = u"m^6/mol^2/s"]
        K_eq = 13.0 / 1000.0,
        [description = "Equilibrium constant in denominator", unit = u"m^3/mol"]
    end

    @variables begin
        H_plus(t), [description = "Hydrogen ion concentration", unit = u"mol/m^3"]
        H2O2_aq(t), [description = "Aqueous H2O2 concentration", unit = u"mol/m^3"]
        HSO3_minus(t), [description = "Bisulfite concentration [HSO3-]", unit = u"mol/m^3"]
        R_H2O2(t), [description = "S(IV) oxidation rate by H2O2", unit = u"mol/m^3/s"]
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

  - H_plus: Hydrogen ion concentration (mol/m^3)
  - Fe_III: Fe(III) concentration (mol/m^3)
  - S_IV_total: Total S(IV) concentration (mol/m^3)
  - R_Fe: Oxidation rate (mol/m^3/s)
"""
@component function SulfateFormationFe(; name = :SulfateFe)
    @constants begin
        k_low_pH = 6.0, [description = "Rate constant for pH 0-3.6", unit = u"s^-1"]
    end

    @variables begin
        H_plus(t), [description = "Hydrogen ion concentration", unit = u"mol/m^3"]
        Fe_III(t), [description = "Fe(III) concentration", unit = u"mol/m^3"]
        S_IV_total(t), [description = "Total S(IV) concentration", unit = u"mol/m^3"]
        R_Fe(t),
        [description = "S(IV) oxidation rate by Fe(III)-catalyzed O2", unit = u"mol/m^3/s"]
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

  - Mn_II: Mn(II) concentration (mol/m^3)
  - S_IV_total: Total S(IV) concentration (mol/m^3)
  - R_Mn: Oxidation rate (mol/m^3/s)
"""
@component function SulfateFormationMn(; name = :SulfateMn)
    @constants begin
        k0 = 1000.0 / 1000.0,
        [description = "Rate constant for low S(IV)", unit = u"m^3/mol/s"]
    end

    @variables begin
        Mn_II(t), [description = "Mn(II) concentration", unit = u"mol/m^3"]
        S_IV_total(t), [description = "Total S(IV) concentration", unit = u"mol/m^3"]
        R_Mn(t),
        [description = "S(IV) oxidation rate by Mn(II)-catalyzed O2", unit = u"mol/m^3/s"]
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

  - Mn_II: Mn(II) concentration (mol/m^3)
  - Fe_III: Fe(III) concentration (mol/m^3)
  - S_IV_total: Total S(IV) concentration (mol/m^3)
  - R_Mn_term: Mn-only contribution (mol/m^3/s)
  - R_Fe_term: Fe-only contribution (mol/m^3/s)
  - R_synergy_term: Synergistic contribution (mol/m^3/s)
  - R_FeMn: Total oxidation rate (mol/m^3/s)
"""
@component function SulfateFormationFeMn(; name = :SulfateFeMn)
    @constants begin
        k_Mn = 750.0 / 1000.0,
        [description = "Rate constant for Mn-only term", unit = u"m^3/mol/s"]
        k_Fe = 2600.0 / 1000.0,
        [description = "Rate constant for Fe-only term", unit = u"m^3/mol/s"]
        k_FeMn = 1.0e10 / 1e6,
        [description = "Rate constant for synergistic term", unit = u"m^6/mol^2/s"]
    end

    @variables begin
        Mn_II(t), [description = "Mn(II) concentration", unit = u"mol/m^3"]
        Fe_III(t), [description = "Fe(III) concentration", unit = u"mol/m^3"]
        S_IV_total(t), [description = "Total S(IV) concentration", unit = u"mol/m^3"]
        R_Mn_term(t), [description = "Mn-only contribution to rate", unit = u"mol/m^3/s"]
        R_Fe_term(t), [description = "Fe-only contribution to rate", unit = u"mol/m^3/s"]
        R_synergy_term(t),
        [description = "Synergistic contribution to rate", unit = u"mol/m^3/s"]
        R_FeMn(t),
        [
            description = "Total S(IV) oxidation rate with Fe/Mn synergism", unit = u"mol/m^3/s"]
    end

    eqs = [
        # Eq 7.102: Three-term rate expression
        R_Mn_term ~ k_Mn * Mn_II * S_IV_total,
        R_Fe_term ~ k_Fe * Fe_III * S_IV_total,
        R_synergy_term ~ k_FeMn * Mn_II * Fe_III * S_IV_total,

        # Total rate
        R_FeMn ~ R_Mn_term + R_Fe_term + R_synergy_term
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

  - R_total: Total S(IV) oxidation rate (mol/m^3/s)
  - S_VI_production: Sulfate production rate (mol/m^3/s)

Also calculates rate conversion factors (Section 7.4):

  - R_ppb_hr: Rate in ppb hr^-1 (Eq 7.75)
  - R_percent_hr: Rate in % hr^-1 (Eq 7.76)
  - tau_SO2: SO2 lifetime (s) (Eq 7.78)    # Create subsystems
"""
@component function SulfateFormation(; name = :SulfateForm)
    # Create subsystems
    o3_ox = SulfateFormationO3(; name = :o3_ox)
    h2o2_ox = SulfateFormationH2O2(; name = :h2o2_ox)
    fe_ox = SulfateFormationFe(; name = :fe_ox)
    mn_ox = SulfateFormationMn(; name = :mn_ox)
    femn_ox = SulfateFormationFeMn(; name = :femn_ox)

    @constants begin
        R_gas = 8.31446,
        [description = "Gas constant for rate conversion", unit = u"J/mol/K"]
        p_total = 101325.0, [description = "Standard atmospheric pressure", unit = u"Pa"]
        # Inverse water density: 1/(10^6 g/m^3) = 10^-6 m^3/g
        rho_w_inv = 1e-6,
        [description = "Inverse water density for LWC conversion", unit = u"m^3/g"]
        # Time conversion: seconds to hours
        s_per_hr = 3600.0, [description = "Seconds per hour (dimensionless)", unit = u"s"]
        # ppb scaling factor (dimensionless)
        ppb_factor = 1e9,
        [description = "Parts per billion scaling factor (dimensionless)", unit = u"1"]
        # percent scaling factor (dimensionless)
        pct_factor = 100.0,
        [description = "Percent scaling factor (dimensionless)", unit = u"1"]
    end

    @parameters begin
        L, [description = "Liquid water content", unit = u"g/m^3"]
        xi_SO2, [description = "SO2 mixing ratio (dimensionless)", unit = u"1"]
    end

    @variables begin
        T(t), [description = "Temperature", unit = u"K"]
        H_plus(t), [description = "Hydrogen ion concentration", unit = u"mol/m^3"]
        SO2_aq(t), [description = "Aqueous SO2 concentration [SO2.H2O]", unit = u"mol/m^3"]
        HSO3_minus(t), [description = "Bisulfite concentration [HSO3-]", unit = u"mol/m^3"]
        SO3_2minus(t), [description = "Sulfite concentration [SO3^2-]", unit = u"mol/m^3"]
        S_IV_total(t), [description = "Total S(IV) concentration", unit = u"mol/m^3"]
        O3_aq(t), [description = "Aqueous O3 concentration", unit = u"mol/m^3"]
        H2O2_aq(t), [description = "Aqueous H2O2 concentration", unit = u"mol/m^3"]
        Fe_III(t), [description = "Fe(III) concentration", unit = u"mol/m^3"]
        Mn_II(t), [description = "Mn(II) concentration", unit = u"mol/m^3"]
        R_total(t), [description = "Total S(IV) oxidation rate", unit = u"mol/m^3/s"]
        S_VI_production(t), [description = "Sulfate production rate", unit = u"mol/m^3/s"]
        w_L(t), [description = "Liquid water volume ratio (dimensionless)", unit = u"1"]
        R_ppb_hr(t), [description = "Rate in ppb/hr (Eq 7.75) (dimensionless)", unit = u"1"]
        R_percent_hr(t),
        [description = "Rate in percent/hr (Eq 7.76) (dimensionless)", unit = u"1"]
        tau_SO2(t), [description = "SO2 lifetime (Eq 7.78)", unit = u"s"]
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

        # Liquid water volume ratio: w_L = L / rho_water (dimensionless)
        w_L ~ rho_w_inv * L,

        # Eq 7.75: Rate conversion to ppb/hr
        # R''_a (ppb h^-1) = 10^9 * w_L * R * T * R_a / p_total * 3600
        # Units: (1) * (Pa·m³/(mol·K)) * K * (mol/(m³·s)) / Pa * s = (1) ✓
        R_ppb_hr ~ ppb_factor * s_per_hr * w_L * R_gas * T * R_total / p_total,

        # Eq 7.76: Rate conversion to %/hr
        # R_bar_a (%h^-1) = 100 * w_L * R * T * R_a / (p_total * xi_SO2) * 3600
        R_percent_hr ~
        pct_factor * s_per_hr * R_total * w_L * R_gas * T / (p_total * xi_SO2),

        # Eq 7.78: SO2 lifetime
        # tau_SO2 (s) = xi_SO2 * p_total / (w_L * R * T * R_a)
        tau_SO2 ~ xi_SO2 * p_total / (w_L * R_gas * T * R_total)
    ]

    return System(eqs, t;
        systems = [o3_ox, h2o2_ox, fe_ox, mn_ox, femn_ox],
        name)
end

# =============================================================================
# Rate Conversion Utilities (Section 7.4)
# =============================================================================

"""
    rate_to_ppb_hr(R_a, L, T; p_total=101325.0)

Convert aqueous-phase rate to atmospheric mixing ratio rate (Eq 7.75).

Arguments:

  - R_a: Aqueous reaction rate (mol/L/s)
  - L: Liquid water content (g/m³)
  - T: Temperature (K)
  - p_total: Total pressure (Pa), default 101325 (1 atm)

Returns:

  - Rate in ppb/hr
"""
function rate_to_ppb_hr(R_a, L, T; p_total = 101325.0)
    R = 8314.46  # Pa L mol^-1 K^-1
    return 3.6e6 * L * R * T * R_a / p_total
end

"""
    rate_to_percent_hr(R_a, L, T, xi_SO2; p_total=101325.0)

Convert aqueous-phase rate to percent conversion rate (Eq 7.76).

Arguments:

  - R_a: Aqueous reaction rate (mol/L/s)
  - L: Liquid water content (g/m³)
  - T: Temperature (K)
  - xi_SO2: SO2 mixing ratio (dimensionless)
  - p_total: Total pressure (Pa), default 101325 (1 atm)

Returns:

  - Rate in %/hr
"""
function rate_to_percent_hr(R_a, L, T, xi_SO2; p_total = 101325.0)
    R = 8314.46  # Pa L mol^-1 K^-1
    return 3.6e8 * R_a * L * R * T / (p_total * xi_SO2)
end

"""
    so2_lifetime(R_a, L, T, xi_SO2; p_total=101325.0)

Calculate SO2 lifetime in seconds (Eq 7.78).

Arguments:

  - R_a: Aqueous reaction rate (mol/L/s)
  - L: Liquid water content (g/m³)
  - T: Temperature (K)
  - xi_SO2: SO2 mixing ratio (dimensionless)
  - p_total: Total pressure (Pa), default 101325 (1 atm)

Returns:

  - SO2 lifetime in seconds
"""
function so2_lifetime(R_a, L, T, xi_SO2; p_total = 101325.0)
    R = 8314.46  # Pa L mol^-1 K^-1
    return xi_SO2 * p_total / (1e3 * R_a * L * R * T)
end
