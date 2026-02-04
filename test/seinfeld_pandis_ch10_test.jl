@testitem "KelvinEffect Structure" tags=[:seinfeld_pandis_ch10] begin
    using Aerosol
    using Aerosol.SeinfeldPandisCh10
    using ModelingToolkit

    # Test that the component builds correctly
    sys = KelvinEffect()
    @test sys isa ModelingToolkit.AbstractSystem

    # Check that expected variables exist
    var_names = [Symbol(ModelingToolkit.getname(v)) for v in unknowns(sys)]
    @test :ln_S in var_names
    @test :S in var_names

    # Check that expected parameters exist
    param_names = [Symbol(ModelingToolkit.getname(p)) for p in parameters(sys)]
    @test :T in param_names
    @test :R_p in param_names
    @test :sigma in param_names || :σ in param_names

    # Verify equation count: 2 equations (ln_S and S)
    @test length(equations(sys)) == 2
end

@testitem "Kelvin Effect Calculation" tags=[:seinfeld_pandis_ch10] begin
    using Aerosol.SeinfeldPandisCh10

    # Eq. 10.86: S = exp(2σM / (R T ρ_l R_p))
    # For a 50 nm water droplet at 298 K:
    # S = exp(2 * 0.072 * 0.018 / (8.314 * 298 * 1000 * 50e-9))
    R_gas = 8.314
    exponent = 2 * 0.072 * 0.018 / (R_gas * 298.0 * 1000.0 * 50e-9)
    S_expected = exp(exponent)

    S = kelvin_saturation_ratio(298.0, 50e-9, 0.072, 0.018, 1000.0)
    @test S ≈ S_expected rtol=1e-10

    # Test limiting cases
    # Large particles: Kelvin effect should be negligible (S → 1)
    S_large = kelvin_saturation_ratio(298.0, 1e-6, 0.072, 0.018, 1000.0)
    @test S_large ≈ 1.0 atol=0.002

    S_very_large = kelvin_saturation_ratio(298.0, 1e-3, 0.072, 0.018, 1000.0)
    @test S_very_large ≈ 1.0 atol=1e-5

    # Small particles: Kelvin effect should be significant
    S_small = kelvin_saturation_ratio(298.0, 10e-9, 0.072, 0.018, 1000.0)
    @test S_small > 1.1

    # Temperature dependence: higher T → lower Kelvin effect
    S_cold = kelvin_saturation_ratio(250.0, 50e-9, 0.072, 0.018, 1000.0)
    S_warm = kelvin_saturation_ratio(350.0, 50e-9, 0.072, 0.018, 1000.0)
    @test S_cold > S > S_warm

    # Surface tension dependence: higher σ → larger Kelvin effect
    S_low_sigma = kelvin_saturation_ratio(298.0, 50e-9, 0.030, 0.018, 1000.0)
    S_high_sigma = kelvin_saturation_ratio(298.0, 50e-9, 0.072, 0.018, 1000.0)
    @test S_high_sigma > S_low_sigma

    # Monotonicity: S always >= 1 for physical parameters
    for R_p in [5e-9, 10e-9, 50e-9, 100e-9, 1e-6]
        @test kelvin_saturation_ratio(298.0, R_p, 0.072, 0.018, 1000.0) >= 1.0
    end

    # Table 10.6: Verify with different compounds
    # Benzene: M=0.07811 kg/mol, ρ=879 kg/m³, σ=0.02821 N/m
    S_benzene = kelvin_saturation_ratio(298.0, 50e-9, 0.02821, 0.07811, 879.0)
    @test S_benzene > 1.0
end

@testitem "DRH Temperature Dependence" tags=[:seinfeld_pandis_ch10] begin
    using Aerosol
    using Aerosol.SeinfeldPandisCh10
    using ModelingToolkit

    # Test the DRH component builds
    sys = DRHTemperature(salt = :NH4NO3)
    @test sys isa ModelingToolkit.AbstractSystem

    # Test structural: 1 equation for DRH
    @test length(equations(sys)) == 1

    # Test functional form at reference temperature
    # NH4NO3 DRH at 298 K should be ~0.618
    drh_298 = drh_temperature(298.0, :NH4NO3)
    @test drh_298 ≈ 0.618 rtol=0.05

    # DRH decreases with increasing temperature for salts with positive ΔH_s
    # (NH4NO3 has ΔH_s = 16.27 kJ/mol > 0)
    drh_280 = drh_temperature(280.0, :NH4NO3)
    drh_310 = drh_temperature(310.0, :NH4NO3)
    @test drh_280 > drh_298 > drh_310

    # NaCl has ΔH_s = 1.88 kJ/mol > 0, so DRH should also decrease with T
    drh_nacl_280 = drh_temperature(280.0, :NaCl)
    drh_nacl_310 = drh_temperature(310.0, :NaCl)
    @test drh_nacl_280 > drh_nacl_310

    # Na2SO4 has ΔH_s = -9.76 kJ/mol < 0, so DRH should increase with T
    drh_na2so4_280 = drh_temperature(280.0, :Na2SO4)
    drh_na2so4_310 = drh_temperature(310.0, :Na2SO4)
    @test drh_na2so4_280 < drh_na2so4_310

    # Test empirical NH4NO3 DRH formula (Eq. 10.88)
    drh_emp = nh4no3_drh(298.0)
    @test drh_emp ≈ 0.618 rtol=0.05

    # DRH must be between 0 and 1 for physical temperatures
    for T in [250.0, 270.0, 298.0, 310.0, 330.0]
        drh_val = drh_temperature(T, :NH4NO3)
        @test 0.0 < drh_val < 1.0
    end

    # Test all salts with data build correctly
    for salt in [:NH42SO4, :NaCl, :NaNO3, :KCl, :Na2SO4, :NH4NO3]
        sys_salt = DRHTemperature(salt = salt)
        @test sys_salt isa ModelingToolkit.AbstractSystem
    end
end

@testitem "DRH Data Tables" tags=[:seinfeld_pandis_ch10] begin
    using Aerosol.SeinfeldPandisCh10

    # Verify all DRH values at 298 K match Table 10.1
    @test SeinfeldPandisCh10.DRH_298[:KCl] ≈ 0.842 rtol=0.01
    @test SeinfeldPandisCh10.DRH_298[:Na2SO4] ≈ 0.842 rtol=0.01
    @test SeinfeldPandisCh10.DRH_298[:NH4Cl] ≈ 0.800 rtol=0.01
    @test SeinfeldPandisCh10.DRH_298[:NH42SO4] ≈ 0.799 rtol=0.01
    @test SeinfeldPandisCh10.DRH_298[:NaCl] ≈ 0.753 rtol=0.01
    @test SeinfeldPandisCh10.DRH_298[:NaNO3] ≈ 0.743 rtol=0.01
    @test SeinfeldPandisCh10.DRH_298[:NH43HSO42] ≈ 0.690 rtol=0.01
    @test SeinfeldPandisCh10.DRH_298[:NH4NO3] ≈ 0.618 rtol=0.01
    @test SeinfeldPandisCh10.DRH_298[:NaHSO4] ≈ 0.520 rtol=0.01
    @test SeinfeldPandisCh10.DRH_298[:NH4HSO4] ≈ 0.400 rtol=0.01

    # Verify enthalpy of solution values match Table 10.3
    @test SeinfeldPandisCh10.DELTA_HS_298[:NH42SO4] ≈ 6.32 rtol=0.01
    @test SeinfeldPandisCh10.DELTA_HS_298[:Na2SO4] ≈ -9.76 rtol=0.01
    @test SeinfeldPandisCh10.DELTA_HS_298[:NaNO3] ≈ 13.24 rtol=0.01
    @test SeinfeldPandisCh10.DELTA_HS_298[:NH4NO3] ≈ 16.27 rtol=0.01
    @test SeinfeldPandisCh10.DELTA_HS_298[:KCl] ≈ 15.34 rtol=0.01
    @test SeinfeldPandisCh10.DELTA_HS_298[:NaCl] ≈ 1.88 rtol=0.01

    # Verify solubility parameters at 298 K match Table 10.2
    for (salt, expected_n) in [(:NH42SO4, 0.104), (:Na2SO4, 0.065), (:NaNO3, 0.194),
        (:NH4NO3, 0.475), (:KCl, 0.086), (:NaCl, 0.111)]
        A, B, C = SeinfeldPandisCh10.SOLUBILITY_PARAMS[salt]
        n_298 = A + B * 298.0 + C * 298.0^2
        @test n_298 ≈ expected_n rtol=0.05
    end
end

@testitem "ZSR Water Content" tags=[:seinfeld_pandis_ch10] begin
    using Aerosol
    using Aerosol.SeinfeldPandisCh10
    using ModelingToolkit

    # Test the ZSR component builds
    sys = ZSRWaterContent(n_species = 2)
    @test sys isa ModelingToolkit.AbstractSystem

    # Verify equation count: 2 equations (α_w and W)
    @test length(equations(sys)) == 2

    # Test functional form (Eq. 10.98)
    concentrations = Dict(:NH42SO4 => 1e-6, :NH4NO3 => 1e-6)
    W = zsr_water_content(0.8, concentrations)
    @test W > 0

    # Higher RH → more water (binary molality decreases at higher a_w)
    W_high = zsr_water_content(0.9, concentrations)
    W_low = zsr_water_content(0.6, concentrations)
    @test W_high > W > W_low

    # More salt → more water (linear in concentration by ZSR)
    W_more = zsr_water_content(0.8, Dict(:NH42SO4 => 2e-6, :NH4NO3 => 2e-6))
    @test W_more ≈ 2 * W rtol=0.01  # ZSR is linear in concentration

    # Zero concentration → zero water
    W_zero = zsr_water_content(0.8, Dict(:NH42SO4 => 0.0, :NH4NO3 => 0.0))
    @test W_zero ≈ 0.0 atol=1e-20

    # Single-component ZSR should equal C/m0(a_w)
    W_single = zsr_water_content(0.8, Dict(:NaCl => 1e-6))
    m0_nacl = SeinfeldPandisCh10.binary_molality_nacl(0.8)
    @test W_single ≈ 1e-6 / m0_nacl rtol=1e-10

    # Test additivity: W(A+B) = W(A) + W(B) (ZSR mixing rule)
    W_a = zsr_water_content(0.8, Dict(:NH42SO4 => 1e-6))
    W_b = zsr_water_content(0.8, Dict(:NH4NO3 => 1e-6))
    W_ab = zsr_water_content(0.8, Dict(:NH42SO4 => 1e-6, :NH4NO3 => 1e-6))
    @test W_ab ≈ W_a + W_b rtol=1e-10

    # Binary molality functions should be positive and monotonically decreasing
    for a_w in [0.5, 0.6, 0.7, 0.8, 0.9]
        @test SeinfeldPandisCh10.binary_molality_nh42so4(a_w) > 0
        @test SeinfeldPandisCh10.binary_molality_nh4no3(a_w) > 0
        @test SeinfeldPandisCh10.binary_molality_nacl(a_w) > 0
    end
    # Higher a_w → lower molality (more dilute solution)
    for func in [SeinfeldPandisCh10.binary_molality_nh42so4,
        SeinfeldPandisCh10.binary_molality_nh4no3,
        SeinfeldPandisCh10.binary_molality_nacl]
        @test func(0.6) > func(0.8) > func(0.95)
    end
end

@testitem "NH4NO3 Equilibrium Constants" tags=[:seinfeld_pandis_ch10] begin
    using Aerosol
    using Aerosol.SeinfeldPandisCh10
    using ModelingToolkit

    # Test the component builds
    sys = NH4NO3Equilibrium()
    @test sys isa ModelingToolkit.AbstractSystem

    # Verify equation count: 3 equations (ln_Kp, DRH_out, is_aqueous)
    @test length(equations(sys)) == 3

    # Test Kp calculation (Eq. 10.91)
    # At 298 K: ln_Kp = 84.6 - 24220/298 - 6.1*ln(298/298) = 84.6 - 81.28 - 0 = 3.32
    ln_Kp_expected = 84.6 - 24220.0/298.0 - 6.1 * log(298.0/298.0)
    @test ln_Kp_expected ≈ 3.32 atol=0.01
    Kp_298 = nh4no3_Kp(298.0)
    @test Kp_298 ≈ exp(ln_Kp_expected) rtol=1e-10

    # Temperature dependence: Kp increases strongly with temperature
    Kp_280 = nh4no3_Kp(280.0)
    Kp_310 = nh4no3_Kp(310.0)
    @test Kp_310 > Kp_298 > Kp_280

    # Kp should be positive for all temperatures
    for T in [250.0, 270.0, 298.0, 310.0, 330.0]
        @test nh4no3_Kp(T) > 0
    end

    # Test K_AN calculation (Eq. 10.97)
    K_AN_298 = nh4no3_K_AN(298.0)
    @test K_AN_298 ≈ 4.0e17 rtol=1e-10  # Exact at reference temperature

    # K_AN temperature dependence
    K_AN_280 = nh4no3_K_AN(280.0)
    K_AN_310 = nh4no3_K_AN(310.0)
    @test K_AN_280 > K_AN_298 > K_AN_310

    # Test ionic strength fraction (Eq. 10.100)
    # Y = C_NH4NO3 / (C_NH4NO3 + 3*C_NH42SO4)
    # Pure NH4NO3: Y = 1/(1+0) = 1
    @test ionic_strength_fraction(1.0, 0.0) ≈ 1.0 atol=1e-10

    # Equal concentrations: Y = 1/(1+3) = 0.25
    @test ionic_strength_fraction(1.0, 1.0) ≈ 0.25 rtol=0.01

    # Dominant sulfate: Y → 0
    @test ionic_strength_fraction(0.01, 10.0) < 0.001

    # Y is always between 0 and 1
    for (c1, c2) in [(0.1, 0.5), (1.0, 0.1), (5.0, 5.0)]
        Y = ionic_strength_fraction(c1, c2)
        @test 0.0 <= Y <= 1.0
    end

    # Test NH4NO3 DRH empirical formula (Eq. 10.88)
    # ln(DRH%) = 723.7/T + 1.6954
    # At 298 K: ln(DRH%) = 723.7/298 + 1.6954 = 2.428 + 1.6954 = 4.124
    # DRH% = exp(4.124) = 61.8, DRH = 0.618
    drh_298 = nh4no3_drh(298.0)
    @test drh_298 ≈ 0.618 rtol=0.01
end

@testitem "Table 10.7 Equilibrium Constants" tags=[:seinfeld_pandis_ch10] begin
    using Aerosol.SeinfeldPandisCh10

    # Verify all 13 equilibrium constants at 298 K match Table 10.7
    @test get_equilibrium_constant(:NaCl_HNO3, 298.0) ≈ 3.96 rtol=1e-6
    @test get_equilibrium_constant(:HSO4_dissoc, 298.0) ≈ 1.01e-2 rtol=1e-6
    @test get_equilibrium_constant(:NH3_HNO3_aq, 298.0) ≈ 4.0e17 rtol=1e-6
    @test get_equilibrium_constant(:HCl_dissoc, 298.0) ≈ 2.03e6 rtol=1e-6
    @test get_equilibrium_constant(:NH3_HCl_aq, 298.0) ≈ 2.12e17 rtol=1e-6
    @test get_equilibrium_constant(:Na2SO4_dissoc, 298.0) ≈ 0.48 rtol=1e-6
    @test get_equilibrium_constant(:NH42SO4_dissoc, 298.0) ≈ 1.425 rtol=1e-6
    @test get_equilibrium_constant(:HNO3_dissoc, 298.0) ≈ 3.638e6 rtol=1e-6
    @test get_equilibrium_constant(:NH4Cl_dissoc, 298.0) ≈ 1.039e-16 rtol=1e-6
    @test get_equilibrium_constant(:NH4NO3_solid, 298.0) ≈ 3.35e16 rtol=1e-6
    @test get_equilibrium_constant(:NaCl_dissoc, 298.0) ≈ 37.74 rtol=1e-6
    @test get_equilibrium_constant(:NaHSO4_dissoc, 298.0) ≈ 2.44e4 rtol=1e-6
    @test get_equilibrium_constant(:NaNO3_dissoc, 298.0) ≈ 11.97 rtol=1e-6

    # Verify temperature dependence formula: K(T) = K(298) * exp{a(298/T-1) + b[1+ln(298/T)-298/T]}
    # At T=298, the exponential factor should be exactly 1
    for rxn in keys(SeinfeldPandisCh10.EQUILIBRIUM_CONSTANTS)
        K_ref = SeinfeldPandisCh10.EQUILIBRIUM_CONSTANTS[rxn][1]
        @test get_equilibrium_constant(rxn, 298.0) ≈ K_ref rtol=1e-10
    end

    # Temperature sensitivity: verify K changes with temperature
    for rxn in keys(SeinfeldPandisCh10.EQUILIBRIUM_CONSTANTS)
        K_280 = get_equilibrium_constant(rxn, 280.0)
        K_298 = get_equilibrium_constant(rxn, 298.0)
        K_310 = get_equilibrium_constant(rxn, 310.0)
        # All K values should be positive
        @test K_280 > 0
        @test K_298 > 0
        @test K_310 > 0
        # K should change with temperature (not be constant)
        @test K_280 != K_298
    end

    # Cross-consistency: NH3_HNO3_aq K at 298 should match K_AN from nh4no3_K_AN
    @test get_equilibrium_constant(:NH3_HNO3_aq, 298.0) ≈ nh4no3_K_AN(298.0) rtol=0.01

    # Error handling: unknown reaction should throw
    @test_throws ErrorException get_equilibrium_constant(:nonexistent, 298.0)
end

@testitem "ModelingToolkit Integration" tags=[:seinfeld_pandis_ch10] begin
    using Aerosol
    using Aerosol.SeinfeldPandisCh10
    using ModelingToolkit
    using ModelingToolkit: mtkcompile
    using OrdinaryDiffEq
    using SciMLBase: ReturnCode

    # Test that Kelvin effect can be compiled and simulated
    sys = KelvinEffect()
    compiled = mtkcompile(sys)

    # Test steady-state solution
    prob = ODEProblem(compiled, [], (0.0, 1.0))
    sol = solve(prob, Tsit5())
    @test sol.retcode == ReturnCode.Success

    # Check that saturation ratio is positive and reasonable
    S_final = sol[compiled.S][end]
    @test S_final > 1.0  # Kelvin effect enhances vapor pressure
    @test S_final < 2.0  # But not by a huge amount for default parameters

    # Test DRH component compilation
    sys_drh = DRHTemperature(salt = :NH4NO3)
    compiled_drh = mtkcompile(sys_drh)
    prob_drh = ODEProblem(compiled_drh, [], (0.0, 1.0))
    sol_drh = solve(prob_drh, Tsit5())
    @test sol_drh.retcode == ReturnCode.Success

    # DRH should be between 0 and 1
    DRH_val = sol_drh[compiled_drh.DRH][end]
    @test 0.0 < DRH_val < 1.0

    # Test ZSR component compilation
    sys_zsr = ZSRWaterContent(n_species = 2)
    compiled_zsr = mtkcompile(sys_zsr)
    prob_zsr = ODEProblem(compiled_zsr,
        Dict(compiled_zsr.C[1] => 1e-6, compiled_zsr.C[2] => 1e-6,
            compiled_zsr.m0[1] => 5.0, compiled_zsr.m0[2] => 5.0),
        (0.0, 1.0))
    sol_zsr = solve(prob_zsr, Tsit5())
    @test sol_zsr.retcode == ReturnCode.Success

    # W should be positive
    W_val = sol_zsr[compiled_zsr.W][end]
    @test W_val > 0

    # Test NH4NO3 component compilation
    sys_nh4 = NH4NO3Equilibrium()
    compiled_nh4 = mtkcompile(sys_nh4)
    prob_nh4 = ODEProblem(compiled_nh4, [], (0.0, 1.0))
    sol_nh4 = solve(prob_nh4, Tsit5())
    @test sol_nh4.retcode == ReturnCode.Success
end
