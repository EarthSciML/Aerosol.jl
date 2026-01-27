"""
Test suite for aqueous chemistry components implementing Seinfeld & Pandis Chapter 7.
"""

@testitem "HenrysLaw_basic" begin
    using ModelingToolkit
    using Aerosol

    sys = HenrysLaw()
    @test sys isa System
    @test length(equations(sys)) == 1
end

@testitem "HenrysLawTemperature" begin
    using ModelingToolkit
    using Aerosol

    sys = HenrysLawTemperature()
    @test sys isa System
    # Should have equations for H_T, C_aq, w_L, f_A, X_aq
    @test length(equations(sys)) >= 5
end

@testitem "Henry_utility_functions" begin
    using Aerosol

    # Test aqueous fraction calculation
    H = 62.0  # NH3 Henry's constant
    R = 0.08205
    T = 298.0
    w_L = 1e-6  # 1 g/m^3 LWC
    X = aqueous_fraction(H, R, T, w_L)
    @test 0 < X < 1

    # Test distribution factor
    f = distribution_factor(H, R, T, w_L)
    @test f > 0
    @test X ≈ f / (1 + f)

    # Test effective Henry's constant
    H_SO2 = 1.23
    K_s1 = 1.3e-2
    K_s2 = 6.6e-8
    H_plus = 1e-4  # pH 4
    H_eff = effective_henrys_constant(H_SO2, K_s1, K_s2, H_plus)
    @test H_eff > H_SO2  # Effective should be larger due to dissociation
end

@testitem "WaterEquilibrium" begin
    using ModelingToolkit
    using Aerosol

    sys = WaterEquilibrium()
    @test sys isa System

    # Check that K_w equation is present
    eqs = equations(sys)
    @test length(eqs) >= 3  # K_w temperature dependence, K_w = H+ * OH-, pH
end

@testitem "CO2Equilibria" begin
    using ModelingToolkit
    using Aerosol

    sys = CO2Equilibria()
    @test sys isa System

    # Should have equations for:
    # - H_CO2, K_c1, K_c2 (temperature dependence)
    # - CO2_aq, HCO3_minus, CO3_2minus (species)
    # - C_total, H_CO2_eff
    eqs = equations(sys)
    @test length(eqs) >= 8
end

@testitem "SO2Equilibria" begin
    using ModelingToolkit
    using Aerosol

    sys = SO2Equilibria()
    @test sys isa System

    # Should have equations for:
    # - H_SO2, K_s1, K_s2 (temperature dependence)
    # - SO2_aq, HSO3_minus, SO3_2minus (species)
    # - S_IV_total, H_SO2_eff
    # - alpha_0, alpha_1, alpha_2 (mole fractions)
    eqs = equations(sys)
    @test length(eqs) >= 11
end

@testitem "NH3Equilibria" begin
    using ModelingToolkit
    using Aerosol

    sys = NH3Equilibria()
    @test sys isa System

    eqs = equations(sys)
    @test length(eqs) >= 7
end

@testitem "HNO3Equilibria" begin
    using ModelingToolkit
    using Aerosol

    sys = HNO3Equilibria()
    @test sys isa System

    eqs = equations(sys)
    @test length(eqs) >= 5
end

@testitem "H2O2Equilibria" begin
    using ModelingToolkit
    using Aerosol

    sys = H2O2Equilibria()
    @test sys isa System

    eqs = equations(sys)
    @test length(eqs) >= 5
end

@testitem "O3Equilibria" begin
    using ModelingToolkit
    using Aerosol

    sys = O3Equilibria()
    @test sys isa System

    eqs = equations(sys)
    @test length(eqs) >= 2
end

@testitem "SulfateFormationO3" begin
    using ModelingToolkit
    using Aerosol

    sys = SulfateFormationO3()
    @test sys isa System

    # Check rate equation
    eqs = equations(sys)
    @test length(eqs) >= 1
end

@testitem "SulfateFormationH2O2" begin
    using ModelingToolkit
    using Aerosol

    sys = SulfateFormationH2O2()
    @test sys isa System

    eqs = equations(sys)
    @test length(eqs) >= 1
end

@testitem "SulfateFormationFe" begin
    using ModelingToolkit
    using Aerosol

    sys = SulfateFormationFe()
    @test sys isa System

    eqs = equations(sys)
    @test length(eqs) >= 1
end

@testitem "SulfateFormationMn" begin
    using ModelingToolkit
    using Aerosol

    sys = SulfateFormationMn()
    @test sys isa System

    eqs = equations(sys)
    @test length(eqs) >= 1
end

@testitem "SulfateFormationFeMn" begin
    using ModelingToolkit
    using Aerosol

    sys = SulfateFormationFeMn()
    @test sys isa System

    eqs = equations(sys)
    @test length(eqs) >= 4  # Mn term, Fe term, synergy term, total
end

@testitem "SulfateFormation_combined" begin
    using ModelingToolkit
    using Aerosol

    sys = SulfateFormation()
    @test sys isa System

    # Should have multiple subsystems
    subsys = ModelingToolkit.get_systems(sys)
    @test length(subsys) == 5  # o3, h2o2, fe, mn, femn
end

@testitem "rate_conversion_functions" begin
    using Aerosol

    # Test rate conversions from Section 7.4
    R_a = 1e-9  # 1 nM/s
    L = 1.0     # 1 g/m^3
    T = 298.0   # K
    xi_SO2 = 1e-9  # 1 ppb

    # Eq 7.75: ppb/hr conversion
    R_ppb = rate_to_ppb_hr(R_a, L, T)
    @test R_ppb > 0

    # Eq 7.76: %/hr conversion
    R_pct = rate_to_percent_hr(R_a, L, T, xi_SO2)
    @test R_pct > 0

    # Eq 7.78: SO2 lifetime
    tau = so2_lifetime(R_a, L, T, xi_SO2)
    @test tau > 0
end

@testitem "physical_constants" begin
    using Aerosol

    # Verify some key constants from the textbook
    @test K_W_298 ≈ 1.0e-14
    @test K_S1_298 ≈ 1.3e-2
    @test K_S2_298 ≈ 6.6e-8
    @test K_C1_298 ≈ 4.3e-7
    @test K_C2_298 ≈ 4.7e-11
    @test K_A1_298 ≈ 1.7e-5

    # Henry's law constants
    @test HENRY_CONSTANTS_298[:SO2] ≈ 1.23
    @test HENRY_CONSTANTS_298[:NH3] ≈ 62.0
    @test HENRY_CONSTANTS_298[:H2O2] ≈ 1.0e5
    @test HENRY_CONSTANTS_298[:HNO3] ≈ 2.1e5
end

@testitem "rate_constants" begin
    using Aerosol

    # Verify S(IV) oxidation rate constants
    @test K0_O3 ≈ 2.4e4
    @test K1_O3 ≈ 3.7e5
    @test K2_O3 ≈ 1.5e9
    @test K_H2O2 ≈ 7.5e7
    @test K_MN_SIV ≈ 750.0
    @test K_FE_SIV ≈ 2600.0
    @test K_FEMN_SIV ≈ 1.0e10
end
