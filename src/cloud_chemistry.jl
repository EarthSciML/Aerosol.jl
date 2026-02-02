"""
    Cloud Chemistry Combined System

Combines all aqueous equilibria and sulfate formation kinetics into a
comprehensive cloud droplet chemistry model from Seinfeld & Pandis Chapter 7.

The main component CloudChemistry includes:
- All gas-liquid equilibria (Henry's law with dissociation)
- pH calculation from electroneutrality (Eq 7.114)
- S(IV) to S(VI) oxidation kinetics
- Dynamic tracking of sulfate production

Key features:
- Temperature-dependent equilibrium constants
- Multiple S(IV) oxidation pathways
- Rate conversion to atmospheric units

All concentrations are in mol/m³, pressures in Pa (SI), temperatures in K.
"""

# =============================================================================
# Cloud Chemistry Combined System
# =============================================================================

"""
    CloudChemistry(; name=:CloudChem)

Comprehensive cloud droplet chemistry system combining all equilibria
and S(IV) oxidation kinetics.

This is the main component for modeling atmospheric aqueous-phase chemistry.
It combines:
1. Water autoionization equilibrium
2. CO2/bicarbonate/carbonate equilibria
3. SO2/bisulfite/sulfite (S(IV)) equilibria
4. NH3/ammonium equilibria
5. HNO3/nitrate equilibria
6. H2O2 equilibria
7. O3 dissolution
8. S(IV) oxidation kinetics (O3, H2O2, Fe, Mn pathways)

The pH is determined from the electroneutrality condition (Eq 7.114):
[H+] + [NH4+] = [OH-] + [HCO3-] + 2[CO3^2-] + [HSO3-] + 2[SO3^2-] + [NO3-] + [Cl-] + 2[SO4^2-]

Parameters:
- L: Liquid water content (g/m³)
- xi_SO2: SO2 mixing ratio for lifetime calculations (dimensionless)

Input Variables (partial pressures and metal concentrations):
- T: Temperature (K)
- p_CO2: CO2 partial pressure (Pa)
- p_SO2: SO2 partial pressure (Pa)
- p_NH3: NH3 partial pressure (Pa)
- p_HNO3: HNO3 partial pressure (Pa)
- p_H2O2: H2O2 partial pressure (Pa)
- p_O3: O3 partial pressure (Pa)
- Fe_III: Fe(III) concentration in droplet (mol/m³)
- Mn_II: Mn(II) concentration in droplet (mol/m³)
- Cl_minus: Chloride concentration (mol/m³)
- SO4_2minus: Sulfate concentration (mol/m³)

Output Variables:
- pH: Droplet pH (dimensionless)
- H_plus: H+ concentration (mol/m³)
- All aqueous species concentrations
- S(IV) oxidation rates
- Sulfate production rate
"""
@component function CloudChemistry(; name=:CloudChem)
    # Create subsystems for equilibria
    water_eq = WaterEquilibrium(; name=:water_eq)
    co2_eq = CO2Equilibria(; name=:co2_eq)
    so2_eq = SO2Equilibria(; name=:so2_eq)
    nh3_eq = NH3Equilibria(; name=:nh3_eq)
    hno3_eq = HNO3Equilibria(; name=:hno3_eq)
    h2o2_eq = H2O2Equilibria(; name=:h2o2_eq)
    o3_eq = O3Equilibria(; name=:o3_eq)

    # Create subsystems for sulfate formation
    sulfate = SulfateFormation(; name=:sulfate)

    @constants begin
        R_gas = 8.31446, [description = "Gas constant", unit = u"J/mol/K"]
        rho_w_inv = 1e-6, [description = "Inverse water density for LWC conversion", unit = u"m^3/g"]
    end

    @parameters begin
        L, [description = "Liquid water content", unit = u"g/m^3"]
        xi_SO2, [description = "SO2 mixing ratio for lifetime calculation (dimensionless)", unit = u"1"]
    end

    @variables begin
        # Primary state variables
        T(t), [description = "Temperature", unit = u"K"]
        H_plus(t), [description = "Hydrogen ion concentration", unit = u"mol/m^3"]
        pH(t), [description = "Droplet pH (dimensionless)", unit = u"1"]

        # Gas-phase partial pressures (inputs)
        p_CO2(t), [description = "Partial pressure of CO2", unit = u"Pa"]
        p_SO2(t), [description = "Partial pressure of SO2", unit = u"Pa"]
        p_NH3(t), [description = "Partial pressure of NH3", unit = u"Pa"]
        p_HNO3(t), [description = "Partial pressure of HNO3", unit = u"Pa"]
        p_H2O2(t), [description = "Partial pressure of H2O2", unit = u"Pa"]
        p_O3(t), [description = "Partial pressure of O3", unit = u"Pa"]

        # Metal catalyst concentrations (inputs)
        Fe_III(t), [description = "Fe(III) concentration in droplet", unit = u"mol/m^3"]
        Mn_II(t), [description = "Mn(II) concentration in droplet", unit = u"mol/m^3"]

        # Background ion concentrations (inputs)
        Cl_minus(t), [description = "Chloride concentration", unit = u"mol/m^3"]
        SO4_2minus(t), [description = "Sulfate concentration", unit = u"mol/m^3"]

        # Derived aqueous concentrations
        OH_minus(t), [description = "Hydroxide concentration", unit = u"mol/m^3"]
        HCO3_minus(t), [description = "Bicarbonate concentration", unit = u"mol/m^3"]
        CO3_2minus(t), [description = "Carbonate concentration", unit = u"mol/m^3"]
        HSO3_minus(t), [description = "Bisulfite concentration", unit = u"mol/m^3"]
        SO3_2minus(t), [description = "Sulfite concentration", unit = u"mol/m^3"]
        NH4_plus(t), [description = "Ammonium concentration", unit = u"mol/m^3"]
        NO3_minus(t), [description = "Nitrate concentration", unit = u"mol/m^3"]
        H2O2_aq(t), [description = "Aqueous H2O2 concentration", unit = u"mol/m^3"]
        O3_aq(t), [description = "Aqueous O3 concentration", unit = u"mol/m^3"]
        SO2_aq(t), [description = "Aqueous SO2 concentration", unit = u"mol/m^3"]

        # Total dissolved species
        S_IV_total(t), [description = "Total S(IV) concentration", unit = u"mol/m^3"]
        C_total(t), [description = "Total dissolved inorganic carbon", unit = u"mol/m^3"]
        NH3_total(t), [description = "Total dissolved ammonia", unit = u"mol/m^3"]
        HNO3_total(t), [description = "Total dissolved nitric acid", unit = u"mol/m^3"]

        # Rate outputs
        R_total(t), [description = "Total S(IV) oxidation rate", unit = u"mol/m^3/s"]
        R_O3(t), [description = "S(IV) oxidation rate by O3", unit = u"mol/m^3/s"]
        R_H2O2(t), [description = "S(IV) oxidation rate by H2O2", unit = u"mol/m^3/s"]
        R_FeMn(t), [description = "S(IV) oxidation rate by Fe/Mn", unit = u"mol/m^3/s"]

        # Charge balance residual (should be zero at equilibrium)
        charge_balance(t), [description = "Electroneutrality residual", unit = u"mol/m^3"]

        # Liquid water mixing ratio
        w_L(t), [description = "Liquid water mixing ratio (dimensionless)", unit = u"1"]
    end

    eqs = [
        # === Temperature coupling ===
        water_eq.T ~ T,
        co2_eq.T ~ T,
        so2_eq.T ~ T,
        nh3_eq.T ~ T,
        hno3_eq.T ~ T,
        h2o2_eq.T ~ T,
        o3_eq.T ~ T,
        sulfate.T ~ T,

        # === H+ coupling (all subsystems share the same H+) ===
        water_eq.H_plus ~ H_plus,
        co2_eq.H_plus ~ H_plus,
        so2_eq.H_plus ~ H_plus,
        nh3_eq.H_plus ~ H_plus,
        hno3_eq.H_plus ~ H_plus,
        h2o2_eq.H_plus ~ H_plus,
        sulfate.H_plus ~ H_plus,

        # === Partial pressure inputs ===
        co2_eq.p_CO2 ~ p_CO2,
        so2_eq.p_SO2 ~ p_SO2,
        nh3_eq.p_NH3 ~ p_NH3,
        hno3_eq.p_HNO3 ~ p_HNO3,
        h2o2_eq.p_H2O2 ~ p_H2O2,
        o3_eq.p_O3 ~ p_O3,

        # === Connect aqueous concentrations from equilibria ===
        OH_minus ~ water_eq.OH_minus,
        HCO3_minus ~ co2_eq.HCO3_minus,
        CO3_2minus ~ co2_eq.CO3_2minus,
        SO2_aq ~ so2_eq.SO2_aq,
        HSO3_minus ~ so2_eq.HSO3_minus,
        SO3_2minus ~ so2_eq.SO3_2minus,
        NH4_plus ~ nh3_eq.NH4_plus,
        NO3_minus ~ hno3_eq.NO3_minus,
        H2O2_aq ~ h2o2_eq.H2O2_aq,
        O3_aq ~ o3_eq.O3_aq,

        # === Total concentrations ===
        S_IV_total ~ so2_eq.S_IV_total,
        C_total ~ co2_eq.C_total,
        NH3_total ~ nh3_eq.NH3_total,
        HNO3_total ~ hno3_eq.HNO3_total,

        # === Connect to sulfate formation subsystem ===
        sulfate.SO2_aq ~ SO2_aq,
        sulfate.HSO3_minus ~ HSO3_minus,
        sulfate.SO3_2minus ~ SO3_2minus,
        sulfate.S_IV_total ~ S_IV_total,
        sulfate.O3_aq ~ O3_aq,
        sulfate.H2O2_aq ~ H2O2_aq,
        sulfate.Fe_III ~ Fe_III,
        sulfate.Mn_II ~ Mn_II,
        sulfate.L ~ L,
        sulfate.xi_SO2 ~ xi_SO2,

        # === Rate outputs from sulfate formation ===
        R_total ~ sulfate.R_total,
        R_O3 ~ sulfate.o3_ox.R_O3,
        R_H2O2 ~ sulfate.h2o2_ox.R_H2O2,
        R_FeMn ~ sulfate.femn_ox.R_FeMn,

        # === pH from water equilibrium ===
        pH ~ water_eq.pH,

        # === Liquid water mixing ratio (Eq 7.1) ===
        w_L ~ rho_w_inv * L,

        # === Electroneutrality (Eq 7.114) ===
        # [H+] + [NH4+] = [OH-] + [HCO3-] + 2[CO3^2-] + [HSO3-] + 2[SO3^2-] + [NO3-] + [Cl-] + 2[SO4^2-]
        # Residual should be zero
        charge_balance ~ (H_plus + NH4_plus) -
                         (OH_minus + HCO3_minus + 2*CO3_2minus +
                          HSO3_minus + 2*SO3_2minus + NO3_minus +
                          Cl_minus + 2*SO4_2minus),
    ]

    return System(eqs, t;
        systems=[water_eq, co2_eq, so2_eq, nh3_eq, hno3_eq, h2o2_eq, o3_eq, sulfate],
        name)
end

# =============================================================================
# Simplified Cloud Chemistry (pH as input)
# =============================================================================

"""
    CloudChemistryFixedpH(; name=:CloudChemFixedpH)

Simplified cloud chemistry system with pH as an input parameter.

This variant takes pH as an input rather than calculating it from
electroneutrality. This is useful for:
- Sensitivity studies varying pH
- Cases where additional ions (not modeled) affect pH
- Quick calculations without solving the full equilibrium

Parameters:
- L: Liquid water content (g m^-3)
- xi_SO2: SO2 mixing ratio for lifetime calculations (dimensionless)

Input Variables:
- T: Temperature (K)
- pH_input: Input pH value (dimensionless)
- All partial pressures and metal concentrations
"""
@component function CloudChemistryFixedpH(; name=:CloudChemFixedpH)
    # Create subsystems for equilibria
    water_eq = WaterEquilibrium(; name=:water_eq)
    co2_eq = CO2Equilibria(; name=:co2_eq)
    so2_eq = SO2Equilibria(; name=:so2_eq)
    nh3_eq = NH3Equilibria(; name=:nh3_eq)
    hno3_eq = HNO3Equilibria(; name=:hno3_eq)
    h2o2_eq = H2O2Equilibria(; name=:h2o2_eq)
    o3_eq = O3Equilibria(; name=:o3_eq)

    # Create subsystems for sulfate formation
    sulfate = SulfateFormation(; name=:sulfate)

    @constants begin
        R_gas = 8.31446, [description = "Gas constant", unit = u"J/mol/K"]
        C_ref = 1000.0, [description = "Reference concentration (1 M = 1000 mol/m³)", unit = u"mol/m^3"]
    end

    @parameters begin
        L, [description = "Liquid water content", unit = u"g/m^3"]
        xi_SO2, [description = "SO2 mixing ratio for lifetime calculation (dimensionless)", unit = u"1"]
    end

    @variables begin
        # Primary state variables
        T(t), [description = "Temperature", unit = u"K"]
        pH_input(t), [description = "Input pH value (dimensionless)", unit = u"1"]
        H_plus(t), [description = "Hydrogen ion concentration", unit = u"mol/m^3"]

        # Gas-phase partial pressures (inputs)
        p_CO2(t), [description = "Partial pressure of CO2", unit = u"Pa"]
        p_SO2(t), [description = "Partial pressure of SO2", unit = u"Pa"]
        p_NH3(t), [description = "Partial pressure of NH3", unit = u"Pa"]
        p_HNO3(t), [description = "Partial pressure of HNO3", unit = u"Pa"]
        p_H2O2(t), [description = "Partial pressure of H2O2", unit = u"Pa"]
        p_O3(t), [description = "Partial pressure of O3", unit = u"Pa"]

        # Metal catalyst concentrations (inputs)
        Fe_III(t), [description = "Fe(III) concentration in droplet", unit = u"mol/m^3"]
        Mn_II(t), [description = "Mn(II) concentration in droplet", unit = u"mol/m^3"]

        # Derived aqueous concentrations
        OH_minus(t), [description = "Hydroxide concentration", unit = u"mol/m^3"]
        HSO3_minus(t), [description = "Bisulfite concentration", unit = u"mol/m^3"]
        SO3_2minus(t), [description = "Sulfite concentration", unit = u"mol/m^3"]
        H2O2_aq(t), [description = "Aqueous H2O2 concentration", unit = u"mol/m^3"]
        O3_aq(t), [description = "Aqueous O3 concentration", unit = u"mol/m^3"]
        SO2_aq(t), [description = "Aqueous SO2 concentration", unit = u"mol/m^3"]

        # Total S(IV)
        S_IV_total(t), [description = "Total S(IV) concentration", unit = u"mol/m^3"]

        # Rate outputs
        R_total(t), [description = "Total S(IV) oxidation rate", unit = u"mol/m^3/s"]
        R_O3(t), [description = "S(IV) oxidation rate by O3", unit = u"mol/m^3/s"]
        R_H2O2(t), [description = "S(IV) oxidation rate by H2O2", unit = u"mol/m^3/s"]
        R_FeMn(t), [description = "S(IV) oxidation rate by Fe/Mn", unit = u"mol/m^3/s"]
    end

    eqs = [
        # === H+ from input pH ===
        H_plus ~ C_ref * 10^(-pH_input),

        # === Temperature coupling ===
        water_eq.T ~ T,
        co2_eq.T ~ T,
        so2_eq.T ~ T,
        nh3_eq.T ~ T,
        hno3_eq.T ~ T,
        h2o2_eq.T ~ T,
        o3_eq.T ~ T,
        sulfate.T ~ T,

        # === H+ coupling ===
        water_eq.H_plus ~ H_plus,
        co2_eq.H_plus ~ H_plus,
        so2_eq.H_plus ~ H_plus,
        nh3_eq.H_plus ~ H_plus,
        hno3_eq.H_plus ~ H_plus,
        h2o2_eq.H_plus ~ H_plus,
        sulfate.H_plus ~ H_plus,

        # === Partial pressure inputs ===
        co2_eq.p_CO2 ~ p_CO2,
        so2_eq.p_SO2 ~ p_SO2,
        nh3_eq.p_NH3 ~ p_NH3,
        hno3_eq.p_HNO3 ~ p_HNO3,
        h2o2_eq.p_H2O2 ~ p_H2O2,
        o3_eq.p_O3 ~ p_O3,

        # === Connect aqueous concentrations ===
        OH_minus ~ water_eq.OH_minus,
        SO2_aq ~ so2_eq.SO2_aq,
        HSO3_minus ~ so2_eq.HSO3_minus,
        SO3_2minus ~ so2_eq.SO3_2minus,
        H2O2_aq ~ h2o2_eq.H2O2_aq,
        O3_aq ~ o3_eq.O3_aq,
        S_IV_total ~ so2_eq.S_IV_total,

        # === Connect to sulfate formation ===
        sulfate.SO2_aq ~ SO2_aq,
        sulfate.HSO3_minus ~ HSO3_minus,
        sulfate.SO3_2minus ~ SO3_2minus,
        sulfate.S_IV_total ~ S_IV_total,
        sulfate.O3_aq ~ O3_aq,
        sulfate.H2O2_aq ~ H2O2_aq,
        sulfate.Fe_III ~ Fe_III,
        sulfate.Mn_II ~ Mn_II,
        sulfate.L ~ L,
        sulfate.xi_SO2 ~ xi_SO2,

        # === Rate outputs ===
        R_total ~ sulfate.R_total,
        R_O3 ~ sulfate.o3_ox.R_O3,
        R_H2O2 ~ sulfate.h2o2_ox.R_H2O2,
        R_FeMn ~ sulfate.femn_ox.R_FeMn,
    ]

    return System(eqs, t;
        systems=[water_eq, co2_eq, so2_eq, nh3_eq, hno3_eq, h2o2_eq, o3_eq, sulfate],
        name)
end

# =============================================================================
# Dynamic Cloud Chemistry ODE System
# =============================================================================

"""
    CloudChemistryODE(; name=:CloudChemODE)

Dynamic cloud chemistry system with time evolution of sulfate.

This system tracks the time evolution of sulfate concentration as S(IV)
is oxidized to S(VI). It can be integrated using OrdinaryDiffEq.

The system includes:
- All equilibria (assumed fast, instantaneous)
- S(IV) oxidation kinetics (rate-limiting)
- Sulfate accumulation

State Variable:
- SO4_2minus: Sulfate concentration (mol/m³), evolves with time

Parameters and other variables same as CloudChemistryFixedpH.
"""
@component function CloudChemistryODE(; name=:CloudChemODE)
    # Create subsystems
    so2_eq = SO2Equilibria(; name=:so2_eq)
    h2o2_eq = H2O2Equilibria(; name=:h2o2_eq)
    o3_eq = O3Equilibria(; name=:o3_eq)
    sulfate = SulfateFormation(; name=:sulfate)

    @constants begin
        R_gas = 8.31446, [description = "Gas constant", unit = u"J/mol/K"]
        C_ref = 1000.0, [description = "Reference concentration (1 M = 1000 mol/m³)", unit = u"mol/m^3"]
    end

    @parameters begin
        L, [description = "Liquid water content", unit = u"g/m^3"]
        xi_SO2, [description = "SO2 mixing ratio (dimensionless)", unit = u"1"]
    end

    @variables begin
        # State variable - sulfate concentration
        SO4_2minus(t), [description = "Sulfate concentration", unit = u"mol/m^3"]

        # Primary inputs
        T(t), [description = "Temperature", unit = u"K"]
        pH_input(t), [description = "Input pH value (dimensionless)", unit = u"1"]
        H_plus(t), [description = "Hydrogen ion concentration", unit = u"mol/m^3"]

        # Gas-phase partial pressures
        p_SO2(t), [description = "Partial pressure of SO2", unit = u"Pa"]
        p_H2O2(t), [description = "Partial pressure of H2O2", unit = u"Pa"]
        p_O3(t), [description = "Partial pressure of O3", unit = u"Pa"]

        # Metal catalysts
        Fe_III(t), [description = "Fe(III) concentration", unit = u"mol/m^3"]
        Mn_II(t), [description = "Mn(II) concentration", unit = u"mol/m^3"]

        # Derived concentrations
        HSO3_minus(t), [description = "Bisulfite concentration", unit = u"mol/m^3"]
        SO3_2minus(t), [description = "Sulfite concentration", unit = u"mol/m^3"]
        SO2_aq(t), [description = "Aqueous SO2 concentration", unit = u"mol/m^3"]
        H2O2_aq(t), [description = "Aqueous H2O2 concentration", unit = u"mol/m^3"]
        O3_aq(t), [description = "Aqueous O3 concentration", unit = u"mol/m^3"]
        S_IV_total(t), [description = "Total S(IV) concentration", unit = u"mol/m^3"]

        # Rate
        R_total(t), [description = "Total S(IV) oxidation rate", unit = u"mol/m^3/s"]
    end

    eqs = [
        # H+ from pH
        H_plus ~ C_ref * 10^(-pH_input),

        # Temperature and H+ coupling
        so2_eq.T ~ T,
        h2o2_eq.T ~ T,
        o3_eq.T ~ T,
        sulfate.T ~ T,

        so2_eq.H_plus ~ H_plus,
        h2o2_eq.H_plus ~ H_plus,
        sulfate.H_plus ~ H_plus,

        # Partial pressures
        so2_eq.p_SO2 ~ p_SO2,
        h2o2_eq.p_H2O2 ~ p_H2O2,
        o3_eq.p_O3 ~ p_O3,

        # Connect aqueous concentrations
        SO2_aq ~ so2_eq.SO2_aq,
        HSO3_minus ~ so2_eq.HSO3_minus,
        SO3_2minus ~ so2_eq.SO3_2minus,
        H2O2_aq ~ h2o2_eq.H2O2_aq,
        O3_aq ~ o3_eq.O3_aq,
        S_IV_total ~ so2_eq.S_IV_total,

        # Connect to sulfate formation
        sulfate.SO2_aq ~ SO2_aq,
        sulfate.HSO3_minus ~ HSO3_minus,
        sulfate.SO3_2minus ~ SO3_2minus,
        sulfate.S_IV_total ~ S_IV_total,
        sulfate.O3_aq ~ O3_aq,
        sulfate.H2O2_aq ~ H2O2_aq,
        sulfate.Fe_III ~ Fe_III,
        sulfate.Mn_II ~ Mn_II,
        sulfate.L ~ L,
        sulfate.xi_SO2 ~ xi_SO2,

        R_total ~ sulfate.R_total,

        # ODE: Sulfate production
        # d[SO4^2-]/dt = R_total (S(IV) oxidation rate)
        D(SO4_2minus) ~ R_total,
    ]

    return System(eqs, t;
        systems=[so2_eq, h2o2_eq, o3_eq, sulfate],
        name)
end
