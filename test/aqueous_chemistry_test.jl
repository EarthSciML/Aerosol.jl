"""
Test suite for aqueous chemistry components implementing Seinfeld & Pandis Chapter 7.
"""

@testsnippet AqueousSetup begin
    using Test
    using ModelingToolkit
    using Aerosol
end

# =============================================================================
# Structural Tests
# =============================================================================

@testitem "HenrysLaw structure" setup=[AqueousSetup] tags=[:aqueous] begin
    sys=HenrysLaw()
    @test sys isa System
    @test length(equations(sys)) == 1
end

@testitem "HenrysLawTemperature structure" setup=[AqueousSetup] tags=[:aqueous] begin
    sys=HenrysLawTemperature()
    @test sys isa System
    # Should have equations for H_T, C_aq, w_L, f_A, X_aq
    @test length(equations(sys)) >= 5
end

@testitem "WaterEquilibrium structure" setup=[AqueousSetup] tags=[:aqueous] begin
    sys=WaterEquilibrium()
    @test sys isa System
    eqs=equations(sys)
    @test length(eqs) >= 3  # K_w temperature dependence, K_w = H+ * OH-, pH
end

@testitem "CO2Equilibria structure" setup=[AqueousSetup] tags=[:aqueous] begin
    sys=CO2Equilibria()
    @test sys isa System
    eqs=equations(sys)
    @test length(eqs) >= 8
end

@testitem "SO2Equilibria structure" setup=[AqueousSetup] tags=[:aqueous] begin
    sys=SO2Equilibria()
    @test sys isa System
    eqs=equations(sys)
    @test length(eqs) >= 11
end

@testitem "NH3Equilibria structure" setup=[AqueousSetup] tags=[:aqueous] begin
    sys=NH3Equilibria()
    @test sys isa System
    eqs=equations(sys)
    @test length(eqs) >= 7
end

@testitem "HNO3Equilibria structure" setup=[AqueousSetup] tags=[:aqueous] begin
    sys=HNO3Equilibria()
    @test sys isa System
    eqs=equations(sys)
    @test length(eqs) >= 5
end

@testitem "H2O2Equilibria structure" setup=[AqueousSetup] tags=[:aqueous] begin
    sys=H2O2Equilibria()
    @test sys isa System
    eqs=equations(sys)
    @test length(eqs) >= 5
end

@testitem "O3Equilibria structure" setup=[AqueousSetup] tags=[:aqueous] begin
    sys=O3Equilibria()
    @test sys isa System
    eqs=equations(sys)
    @test length(eqs) >= 2
end

@testitem "SulfateFormation structure" setup=[AqueousSetup] tags=[:aqueous] begin
    for F in [SulfateFormationO3, SulfateFormationH2O2, SulfateFormationFe,
        SulfateFormationMn, SulfateFormationFeMn]
        sys=F()
        @test sys isa System
        @test length(equations(sys)) >= 1
    end

    # Combined system
    sys=SulfateFormation()
    @test sys isa System
    subsys=ModelingToolkit.get_systems(sys)
    @test length(subsys) == 5  # o3, h2o2, fe, mn, femn
end

@testitem "CloudChemistry structure" setup=[AqueousSetup] tags=[:aqueous] begin
    for F in [CloudChemistry, CloudChemistryFixedpH, CloudChemistryODE]
        sys=F()
        @test sys isa System
    end
end

# =============================================================================
# Physical Constants Verification
# =============================================================================

@testitem "physical constants" setup=[AqueousSetup] tags=[:aqueous] begin
    # Equilibrium constants (mol/L) — unchanged by unit conversion
    @test K_W_298 ≈ 1.0e-14
    @test K_S1_298 ≈ 1.3e-2
    @test K_S2_298 ≈ 6.6e-8
    @test K_C1_298 ≈ 4.3e-7
    @test K_C2_298 ≈ 4.7e-11
    @test K_A1_298 ≈ 1.7e-5

    # Henry's law constants are now in M/Pa (converted from M/atm by dividing by 101325)
    @test HENRY_CONSTANTS_298[:SO2] ≈ 1.23 / ATM_TO_PA
    @test HENRY_CONSTANTS_298[:NH3] ≈ 62.0 / ATM_TO_PA
    @test HENRY_CONSTANTS_298[:H2O2] ≈ 1.0e5 / ATM_TO_PA
    @test HENRY_CONSTANTS_298[:HNO3] ≈ 2.1e5 / ATM_TO_PA

    # Heats of dissolution are now in J/mol (converted from kcal/mol)
    @test DELTA_H_DISSOLUTION[:SO2] ≈ -6.25 * 4184
    @test DELTA_H_DISSOLUTION[:CO2] ≈ -4.85 * 4184
end

@testitem "rate constants" setup=[AqueousSetup] tags=[:aqueous] begin
    # Rate constants are in mol/L-based units (unchanged)
    @test K0_O3 ≈ 2.4e4
    @test K1_O3 ≈ 3.7e5
    @test K2_O3 ≈ 1.5e9
    @test K_H2O2 ≈ 7.5e7
    @test K_MN_SIV ≈ 750.0
    @test K_FE_SIV ≈ 2600.0
    @test K_FEMN_SIV ≈ 1.0e10
end

# =============================================================================
# Equation Verification Tests
# =============================================================================

@testitem "Henry utility functions" setup=[AqueousSetup] tags=[:aqueous] begin
    # Test aqueous fraction calculation with SI units
    H=62.0/ATM_TO_PA  # NH3 Henry's constant in M/Pa
    R=R_GAS_PA           # Pa L mol^-1 K^-1
    T=298.0
    w_L=1e-6  # 1 g/m^3 LWC
    X=aqueous_fraction(H, R, T, w_L)
    @test 0 < X < 1

    # Test distribution factor
    f=distribution_factor(H, R, T, w_L)
    @test f > 0
    @test X ≈ f / (1 + f)

    # Test effective Henry's constant (uses M/Pa internally)
    H_SO2=1.23/ATM_TO_PA  # M/Pa
    K_s1=1.3e-2
    K_s2=6.6e-8
    H_plus=1e-4  # pH 4
    H_eff=effective_henrys_constant(H_SO2, K_s1, K_s2, H_plus)
    @test H_eff > H_SO2  # Effective should be larger due to dissociation

    # Verify effective Henry's constant formula at pH 4 for SO2
    # H* = H * (1 + K_s1/[H+] + K_s1*K_s2/[H+]^2)
    expected=H_SO2*(1+K_s1/H_plus+K_s1*K_s2/H_plus^2)
    @test H_eff ≈ expected
end

@testitem "van't Hoff temperature dependence" setup=[AqueousSetup] tags=[:aqueous] begin
    # Test that Henry's law constant at 298 K returns the original value
    H_298=HENRY_CONSTANTS_298[:SO2]
    dH=DELTA_H_DISSOLUTION[:SO2]
    @test henrys_constant_at_T(H_298, dH, 298.0) ≈ H_298 rtol=1e-10

    # Test that dissolution is exothermic: H increases at lower T for negative dH
    H_280=henrys_constant_at_T(H_298, dH, 280.0)
    @test H_280 > H_298  # For negative dH, H increases at lower T (more soluble)

    # Same check for CO2
    H_CO2_298=HENRY_CONSTANTS_298[:CO2]
    dH_CO2=DELTA_H_DISSOLUTION[:CO2]
    H_CO2_280=henrys_constant_at_T(H_CO2_298, dH_CO2, 280.0)
    @test H_CO2_280 > H_CO2_298
end

@testitem "rate conversion functions" setup=[AqueousSetup] tags=[:aqueous] begin
    # Test rate conversions from Section 7.4
    R_a=1e-9  # 1 nM/s
    L=1.0     # 1 g/m^3
    T=298.0   # K
    xi_SO2=1e-9  # 1 ppb

    # Eq 7.75: ppb/hr conversion
    R_ppb=rate_to_ppb_hr(R_a, L, T)
    @test R_ppb > 0

    # Eq 7.76: %/hr conversion
    R_pct=rate_to_percent_hr(R_a, L, T, xi_SO2)
    @test R_pct > 0

    # Eq 7.78: SO2 lifetime
    tau=so2_lifetime(R_a, L, T, xi_SO2)
    @test tau > 0

    # Consistency check between rate conversions:
    # R_pct = 3.6e8 * R_a * L * R * T / (p_total * xi_SO2)
    # R_ppb = 3.6e6 * L * R * T * R_a / p_total
    # So R_pct / R_ppb = 100 / xi_SO2
    @test R_pct ≈ 100.0 / xi_SO2 * R_ppb rtol=1e-6
end

# =============================================================================
# Conservation Law Tests
# =============================================================================

@testitem "S(IV) alpha values sum to 1" setup=[AqueousSetup] tags=[:aqueous] begin
    # The mole fractions alpha_0 + alpha_1 + alpha_2 should sum to 1
    # for any pH and temperature
    K_s1=K_S1_298
    K_s2=K_S2_298

    for pH in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        H_plus=10.0^(-pH)
        denom=H_plus^2+K_s1*H_plus+K_s1*K_s2
        alpha_0=H_plus^2/denom
        alpha_1=K_s1*H_plus/denom
        alpha_2=K_s1*K_s2/denom
        @test alpha_0 + alpha_1 + alpha_2 ≈ 1.0 atol=1e-12
    end
end

@testitem "Total S(IV) conservation" setup=[AqueousSetup] tags=[:aqueous] begin
    # S_IV_total = SO2_aq + HSO3_minus + SO3_2minus
    # Verify via direct calculation at known conditions
    H_SO2=1.23/ATM_TO_PA  # M/Pa
    K_s1=K_S1_298
    K_s2=K_S2_298
    p_SO2=1e-8*ATM_TO_PA  # 10 ppb in Pa

    for pH in [2.0, 4.0, 6.0]
        H_plus=10.0^(-pH)
        SO2_aq=H_SO2*p_SO2
        HSO3_m=K_s1*SO2_aq/H_plus
        SO3_2m=K_s2*HSO3_m/H_plus
        S_IV=SO2_aq+HSO3_m+SO3_2m
        # Also via effective Henry's law:
        S_IV_via_Heff=H_SO2*p_SO2*(1+K_s1/H_plus+K_s1*K_s2/H_plus^2)
        @test S_IV ≈ S_IV_via_Heff rtol=1e-10
    end
end

# =============================================================================
# Analytical Solution Tests
# =============================================================================

@testitem "SO2 equilibrium concentrations at 298K" setup=[AqueousSetup] tags=[:aqueous] begin
    # Verify SO2 aqueous concentrations at known conditions (T=298K, pH=4)
    # Using direct calculation from equilibrium constants and Henry's law
    H_SO2=1.23/ATM_TO_PA  # M/Pa
    K_s1=K_S1_298  # 1.3e-2 M
    K_s2=K_S2_298  # 6.6e-8 M
    pH=4.0
    H_plus=10.0^(-pH)  # 1e-4 M
    p_SO2=1e-9*ATM_TO_PA  # 1 ppb in Pa

    # Direct equilibrium calculations
    SO2_aq=H_SO2*p_SO2  # mol/L
    HSO3_m=K_s1*SO2_aq/H_plus
    SO3_2m=K_s2*HSO3_m/H_plus

    # At pH 4, K_s1/H_plus = 130 >> 1, so HSO3- dominates over SO2.H2O
    @test HSO3_m > SO2_aq
    # At pH 4, K_s2/H_plus = 6.6e-4 << 1, so SO3^2- << HSO3-
    @test SO3_2m < HSO3_m

    # Effective Henry's constant
    H_eff=H_SO2*(1+K_s1/H_plus+K_s1*K_s2/H_plus^2)
    @test H_eff > H_SO2
    # At pH 4: enhancement factor ≈ 1 + 130 + 0.086 ≈ 131
    enhancement=1+K_s1/H_plus+K_s1*K_s2/H_plus^2
    @test enhancement ≈ 131.086 rtol=1e-3
end

@testitem "CO2 equilibrium at 298K" setup=[AqueousSetup] tags=[:aqueous] begin
    # CO2 system at 298 K, pH 5.6 (natural rain pH)
    H_CO2=3.4e-2/ATM_TO_PA  # M/Pa
    K_c1=K_C1_298  # 4.3e-7 M
    K_c2=K_C2_298  # 4.7e-11 M
    p_CO2=380e-6*ATM_TO_PA  # 380 ppm in Pa
    H_plus=10.0^(-5.6)

    CO2_aq=H_CO2*p_CO2
    HCO3_m=K_c1*CO2_aq/H_plus
    CO3_2m=K_c2*HCO3_m/H_plus

    # At pH 5.6: K_c1/H_plus = 4.3e-7 / 2.51e-6 ≈ 0.17
    # So CO2.H2O dominates, with some HCO3-
    @test CO2_aq > HCO3_m
    @test HCO3_m > CO3_2m  # Carbonate negligible at this pH
end

# =============================================================================
# Limiting Behavior Tests
# =============================================================================

@testitem "S(IV) limiting behavior at extreme pH" setup=[AqueousSetup] tags=[:aqueous] begin
    K_s1=K_S1_298
    K_s2=K_S2_298

    # Very low pH (pH=1): SO2.H2O dominates (alpha_0 → 1)
    H_plus_low=10.0^(-1.0)
    denom_low=H_plus_low^2+K_s1*H_plus_low+K_s1*K_s2
    alpha_0_low=H_plus_low^2/denom_low
    @test alpha_0_low > 0.85  # SO2.H2O is dominant

    # High pH (pH=8): sulfite (SO3^2-) dominates (alpha_2 → 1)
    H_plus_high=10.0^(-8.0)
    denom_high=H_plus_high^2+K_s1*H_plus_high+K_s1*K_s2
    alpha_2_high=K_s1*K_s2/denom_high
    @test alpha_2_high > 0.5  # SO3^2- is dominant at pH 8

    # Intermediate pH (pH=4): HSO3- dominates
    H_plus_mid=10.0^(-4.0)
    denom_mid=H_plus_mid^2+K_s1*H_plus_mid+K_s1*K_s2
    alpha_1_mid=K_s1*H_plus_mid/denom_mid
    @test alpha_1_mid > 0.9  # HSO3- dominates at pH 4
end

@testitem "Effective Henry's constant limiting behavior" setup=[AqueousSetup] tags=[:aqueous] begin
    # H_eff should always be >= H_intrinsic
    H_SO2=1.23/ATM_TO_PA
    K_s1=K_S1_298
    K_s2=K_S2_298

    for pH in [1.0, 3.0, 5.0, 7.0]
        H_plus=10.0^(-pH)
        H_eff=effective_henrys_constant(H_SO2, K_s1, K_s2, H_plus)
        @test H_eff >= H_SO2

        # Enhancement factor should increase with pH (more dissociation)
        if pH>1.0
            H_plus_lower=10.0^(-(pH-1.0))
            H_eff_lower_pH=effective_henrys_constant(H_SO2, K_s1, K_s2, H_plus_lower)
            @test H_eff > H_eff_lower_pH
        end
    end
end

# =============================================================================
# Qualitative Behavior Tests
# =============================================================================

@testitem "H2O2 pathway pH independence" setup=[AqueousSetup] tags=[:aqueous] begin
    # The H2O2 oxidation rate is relatively pH-independent because
    # R = k * [H+] * [H2O2] * [HSO3-] / (1 + K_eq * [H+])
    # and [HSO3-] ∝ 1/[H+], so the [H+] terms partially cancel.

    H_SO2=1.23/ATM_TO_PA
    K_s1=K_S1_298
    p_SO2=1e-9*ATM_TO_PA  # 1 ppb SO2

    H_H2O2=1.0e5/ATM_TO_PA
    p_H2O2=1e-9*ATM_TO_PA  # 1 ppb H2O2

    k=7.5e7  # M^-2 s^-1
    K_eq=13.0  # M^-1

    rates=Float64[]
    for pH in [3.0, 4.0, 5.0]
        H_plus=10.0^(-pH)
        SO2_aq=H_SO2*p_SO2
        HSO3_m=K_s1*SO2_aq/H_plus
        H2O2_aq=H_H2O2*p_H2O2
        R_H2O2=k*H_plus*H2O2_aq*HSO3_m/(1+K_eq*H_plus)
        push!(rates, R_H2O2)
    end
    # Rate varies by less than 2 orders of magnitude over 2 pH units
    @test maximum(rates) / minimum(rates) < 100
end

@testitem "O3 pathway pH dependence" setup=[AqueousSetup] tags=[:aqueous] begin
    # O3 oxidation rate increases strongly with pH because the SO3^2- pathway
    # has a very large rate constant (k2 = 1.5e9 vs k0 = 2.4e4)

    H_SO2=1.23/ATM_TO_PA
    K_s1=K_S1_298
    K_s2=K_S2_298
    p_SO2=1e-9*ATM_TO_PA
    H_O3=1.1e-2/ATM_TO_PA
    p_O3=50e-9*ATM_TO_PA  # 50 ppb O3

    k0=K0_O3
    k1=K1_O3
    k2=K2_O3

    R_low_pH=let
        H_plus=10.0^(-3.0)
        SO2_aq=H_SO2*p_SO2
        HSO3_m=K_s1*SO2_aq/H_plus
        SO3_2m=K_s2*HSO3_m/H_plus
        O3_aq=H_O3*p_O3
        (k0*SO2_aq+k1*HSO3_m+k2*SO3_2m)*O3_aq
    end

    R_high_pH=let
        H_plus=10.0^(-6.0)
        SO2_aq=H_SO2*p_SO2
        HSO3_m=K_s1*SO2_aq/H_plus
        SO3_2m=K_s2*HSO3_m/H_plus
        O3_aq=H_O3*p_O3
        (k0*SO2_aq+k1*HSO3_m+k2*SO3_2m)*O3_aq
    end

    # O3 rate at pH 6 should be much larger than at pH 3
    @test R_high_pH > R_low_pH * 1e3
end

@testitem "Temperature effect on solubility" setup=[AqueousSetup] tags=[:aqueous] begin
    # For all species with negative dH_dissolution, solubility increases at lower T
    for species in [:CO2, :NH3, :SO2, :H2O2, :O3]
        H_298=HENRY_CONSTANTS_298[species]
        dH=DELTA_H_DISSOLUTION[species]
        @test dH < 0  # All dissolution enthalpies should be negative (exothermic)

        H_280=henrys_constant_at_T(H_298, dH, 280.0)
        H_310=henrys_constant_at_T(H_298, dH, 310.0)
        @test H_280 > H_298 > H_310  # More soluble at lower T
    end
end

@testitem "Fe/Mn synergism" setup=[AqueousSetup] tags=[:aqueous] begin
    # The synergistic rate should exceed the sum of individual Fe and Mn rates
    S_IV=1e-6  # 1 µM S(IV)
    Fe=1e-6    # 1 µM Fe(III)
    Mn=1e-6    # 1 µM Mn(II)

    R_Fe_only=K_FE_SIV*Fe*S_IV
    R_Mn_only=K_MN_SIV*Mn*S_IV
    R_synergy=K_FEMN_SIV*Mn*Fe*S_IV
    R_total=R_Fe_only+R_Mn_only+R_synergy

    # Synergistic term should dominate when both metals are present
    @test R_synergy > R_Fe_only + R_Mn_only
end
