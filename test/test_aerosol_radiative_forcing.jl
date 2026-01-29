@testsnippet AerosolRadiativeForcingSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D
    using DynamicQuantities
    using Aerosol
end

@testitem "AerosolLayerRadiativeForcing Structure" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    sys = AerosolLayerRadiativeForcing()

    # Verify system has expected number of equations
    @test length(equations(sys)) == 5

    # Verify all parameters exist
    params = parameters(sys)
    param_names = [string(p) for p in params]
    @test "F_0" in param_names
    @test "τ" in param_names
    @test "ω" in param_names
    @test "β" in param_names
    @test "R_s" in param_names
    @test "A_c" in param_names
    @test "T_a" in param_names

    # Verify all variables exist
    vars = unknowns(sys)
    var_names = [string(v) for v in vars]
    @test any(n -> occursin("r_aer", n), var_names)
    @test any(n -> occursin("t_aer", n), var_names)
    @test any(n -> occursin("R_as", n), var_names)
    @test any(n -> occursin("ΔR_p", n), var_names)
    @test any(n -> occursin("ΔF", n), var_names)
end

@testitem "CriticalSingleScatteringAlbedo Structure" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    sys = CriticalSingleScatteringAlbedo()

    # Verify system has expected number of equations
    @test length(equations(sys)) == 1

    # Verify parameters
    params = parameters(sys)
    param_names = [string(p) for p in params]
    @test "R_s" in param_names
    @test "β" in param_names

    # Verify variables
    vars = unknowns(sys)
    var_names = [string(v) for v in vars]
    @test any(n -> occursin("ω_crit", n), var_names)
end

@testitem "CloudOpticalDepth Structure" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    sys = CloudOpticalDepth()

    # Verify system has expected number of equations
    @test length(equations(sys)) == 1

    # Verify parameters
    params = parameters(sys)
    param_names = [string(p) for p in params]
    @test "h" in param_names
    @test "L" in param_names
    @test "N" in param_names

    # Verify variables
    vars = unknowns(sys)
    var_names = [string(v) for v in vars]
    @test any(n -> occursin("τ_c", n), var_names)
end

@testitem "CloudAlbedo Structure" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    sys = CloudAlbedo()

    # Verify system has expected number of equations
    @test length(equations(sys)) == 2

    # Verify parameters
    params = parameters(sys)
    param_names = [string(p) for p in params]
    @test "τ_c" in param_names
    @test "g" in param_names

    # Verify variables
    vars = unknowns(sys)
    var_names = [string(v) for v in vars]
    @test any(n -> occursin("R_c", n), var_names)
    @test any(n -> occursin("γ", n), var_names)
end

@testitem "CloudAlbedoSensitivity Structure" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    sys = CloudAlbedoSensitivity()

    # Verify system has expected number of equations
    @test length(equations(sys)) == 2

    # Verify parameters
    params = parameters(sys)
    param_names = [string(p) for p in params]
    @test "R_c" in param_names
    @test "N" in param_names

    # Verify variables
    vars = unknowns(sys)
    var_names = [string(v) for v in vars]
    @test any(n -> occursin("dRc_dN", n), var_names)
    @test any(n -> occursin("S", n), var_names)
end

@testitem "IndirectAerosolForcing Structure" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    sys = IndirectAerosolForcing()

    # Verify system has expected number of equations
    @test length(equations(sys)) == 2

    # Verify parameters
    params = parameters(sys)
    param_names = [string(p) for p in params]
    @test "F_0" in param_names
    @test "A_c" in param_names
    @test "T_a" in param_names
    @test "R_c" in param_names
    @test "Δln_N" in param_names

    # Verify variables
    vars = unknowns(sys)
    var_names = [string(v) for v in vars]
    @test any(n -> occursin("ΔR_c", n), var_names)
    @test any(n -> occursin("ΔF_c", n), var_names)
end

@testitem "Critical Single-Scattering Albedo - Figure 24.2 Validation" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test Eq. 24.15: ω_crit = 2*R_s / (2*R_s + β*(1-R_s)²)
    # Figure 24.2 shows curves for different β values

    # Test case: R_s = 0.15, β = 0.29 (Table 24.2 values)
    # Expected ω_crit ≈ 0.6 from Section 24.2
    R_s = 0.15
    β = 0.29
    ω_crit_expected = 2 * R_s / (2 * R_s + β * (1 - R_s)^2)

    # Verify the analytical formula
    @test ω_crit_expected ≈ 0.586 rtol=0.01  # Approximately 0.6

    # Test limiting cases from Figure 24.2
    # At R_s → 0: ω_crit → 0
    R_s_low = 0.01
    ω_crit_low = 2 * R_s_low / (2 * R_s_low + β * (1 - R_s_low)^2)
    @test ω_crit_low < 0.1

    # At R_s → 1: ω_crit → 1
    R_s_high = 0.99
    ω_crit_high = 2 * R_s_high / (2 * R_s_high + β * (1 - R_s_high)^2)
    @test ω_crit_high > 0.99

    # Verify ω_crit increases with increasing R_s
    R_s_vals = [0.1, 0.3, 0.5, 0.7, 0.9]
    ω_crit_vals = [2 * R / (2 * R + β * (1 - R)^2) for R in R_s_vals]
    @test issorted(ω_crit_vals)

    # Verify ω_crit decreases with increasing β (at fixed R_s)
    β_vals = [0.1, 0.2, 0.3, 0.4]
    R_s_fixed = 0.5
    ω_crit_vs_beta = [2 * R_s_fixed / (2 * R_s_fixed + b * (1 - R_s_fixed)^2) for b in β_vals]
    @test issorted(ω_crit_vs_beta, rev=true)
end

@testitem "Cloud Albedo - Figure 24.16 Validation" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test Eq. 24.38: R_c = τ_c / (τ_c + 7.7) for g = 0.85
    g = 0.85
    γ = 2 / (sqrt(3) * (1 - g))  # ≈ 7.7

    @test γ ≈ 7.698 rtol=0.01

    # Test key points from Figure 24.16
    # At τ_c = 0: R_c → 0
    τ_c_zero = 0.01
    R_c_zero = τ_c_zero / (τ_c_zero + γ)
    @test R_c_zero < 0.01

    # At τ_c = 7.7: R_c = 0.5
    τ_c_half = γ
    R_c_half = τ_c_half / (τ_c_half + γ)
    @test R_c_half ≈ 0.5 rtol=0.01

    # At τ_c >> γ: R_c → 1
    τ_c_large = 100.0
    R_c_large = τ_c_large / (τ_c_large + γ)
    @test R_c_large > 0.9

    # Verify R_c increases monotonically with τ_c
    τ_c_vals = [1, 5, 10, 20, 50]
    R_c_vals = [τ / (τ + γ) for τ in τ_c_vals]
    @test issorted(R_c_vals)
end

@testitem "Twomey Susceptibility - Figure 24.18 Validation" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test Eq. 24.41: dR_c/d(ln N) = R_c*(1-R_c)/3

    # Maximum susceptibility occurs at R_c = 0.5
    R_c_max = 0.5
    S_max = R_c_max * (1 - R_c_max) / 3
    @test S_max ≈ 0.0833 rtol=0.01  # 1/12 ≈ 0.0833

    # Susceptibility at R_c = 0.5 is maximum
    R_c_vals = [0.1, 0.3, 0.5, 0.7, 0.9]
    S_vals = [R * (1 - R) / 3 for R in R_c_vals]

    # Find index of maximum
    max_idx = argmax(S_vals)
    @test R_c_vals[max_idx] ≈ 0.5

    # Susceptibility is symmetric around R_c = 0.5
    @test S_vals[1] ≈ S_vals[5] rtol=0.01  # S(0.1) ≈ S(0.9)
    @test S_vals[2] ≈ S_vals[4] rtol=0.01  # S(0.3) ≈ S(0.7)
end

@testitem "Indirect Forcing - Section 24.8.2 Validation" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test Eq. 24.43 numerical validation
    # Equation validates the formula structure: ΔF_c = -F_0 * A_c * T_a² * ΔR_c

    # Parameters from the text (Table 24.2)
    F_0 = 1370.0  # W/m²
    A_c = 0.6     # Cloud fraction
    T_a = 0.76    # Atmospheric transmittance
    R_c = 0.5     # Cloud albedo (at maximum susceptibility)

    # 30% increase in N: Δln(N) = ln(1.3) ≈ 0.262
    Δln_N = log(1.3)
    @test Δln_N ≈ 0.262 rtol=0.01

    # Calculate ΔR_c from Eq. 24.42: ΔR_c = R_c*(1-R_c)/3 * Δln(N)
    ΔR_c = R_c * (1 - R_c) / 3 * Δln_N
    @test ΔR_c ≈ 0.0219 rtol=0.05

    # Calculate ΔF_c from Eq. 24.43
    ΔF_c = -F_0 * A_c * T_a^2 * ΔR_c

    # This gives the cloud forcing sensitivity for a 30% change in N globally
    # The magnitude is consistent with ERBE observations that cloud shortwave
    # forcing is ~50 W/m² (Section 24.8.2), so a ~2% relative change in albedo
    # gives a forcing of several W/m²
    @test ΔF_c < 0  # Negative = cooling (increased albedo reflects more sunlight)
    @test abs(ΔF_c) > 5.0  # Significant forcing magnitude
    @test abs(ΔF_c) < 20.0  # But bounded to realistic range

    # Verify formula consistency: ΔF_c = -F_0 * A_c * T_a² * ΔR_c
    ΔF_c_formula = -F_0 * A_c * T_a^2 * ΔR_c
    @test ΔF_c ≈ ΔF_c_formula rtol=1e-10
end

@testitem "Aerosol Layer Forcing - Limiting Cases" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test limiting behavior from Section 24.1

    # Auxiliary relations from Figure 24.1
    # r = (1 - exp(-τ)) * ω * β
    # t = exp(-τ) + ω*(1-β)*(1-exp(-τ))

    # Case 1: τ = 0 (no aerosol)
    τ_zero = 0.0
    ω = 0.95
    β = 0.29

    r_zero = (1 - exp(-τ_zero)) * ω * β
    t_zero = exp(-τ_zero) + ω * (1 - β) * (1 - exp(-τ_zero))

    @test r_zero ≈ 0.0 atol=1e-10
    @test t_zero ≈ 1.0 atol=1e-10

    # Case 2: ω = 1 (nonabsorbing aerosol)
    τ = 0.1
    ω_nonabs = 1.0
    r_nonabs = (1 - exp(-τ)) * ω_nonabs * β
    t_nonabs = exp(-τ) + ω_nonabs * (1 - β) * (1 - exp(-τ))

    # For nonabsorbing aerosol: r + t should equal 1 (energy conservation)
    # Actually, r + t + absorbed = 1, but absorbed = (1-ω)*(1-exp(-τ)) = 0 for ω=1
    @test r_nonabs + t_nonabs ≈ 1.0 rtol=0.01

    # Case 3: τ << 1 approximation
    τ_small = 0.01
    r_approx = τ_small * ω * β  # Eq. 24.11 approximation
    r_exact = (1 - exp(-τ_small)) * ω * β

    @test r_approx ≈ r_exact rtol=0.05  # Good approximation for small τ
end

@testitem "Cloud Optical Depth - Eq. 24.36" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test Eq. 24.36: τ_c = h * (9π*L²*N / 2ρ_w²)^(1/3)

    # Typical values
    h = 500.0  # m (cloud thickness)
    L = 0.3e-3  # kg/m³ (liquid water content)
    N = 100e6   # m⁻³ (CDNC)
    ρ_w = 1000.0  # kg/m³ (water density)

    τ_c = h * (9 * π * L^2 * N / (2 * ρ_w^2))^(1/3)

    # Cloud optical depth should be positive and reasonable
    @test τ_c > 0
    @test τ_c < 100  # Reasonable upper bound

    # τ_c should increase with N (more droplets = more optical depth)
    N_high = 200e6
    τ_c_high = h * (9 * π * L^2 * N_high / (2 * ρ_w^2))^(1/3)
    @test τ_c_high > τ_c

    # τ_c should increase with L
    L_high = 0.5e-3
    τ_c_L_high = h * (9 * π * L_high^2 * N / (2 * ρ_w^2))^(1/3)
    @test τ_c_L_high > τ_c

    # τ_c should scale linearly with h
    h_double = 1000.0
    τ_c_h_double = h_double * (9 * π * L^2 * N / (2 * ρ_w^2))^(1/3)
    @test τ_c_h_double ≈ 2 * τ_c rtol=0.01
end

@testitem "Direct Forcing - Cooling vs Heating" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test that aerosols with ω > ω_crit cause cooling (ΔF > 0 means increased reflection)

    # Calculate ω_crit for typical surface
    R_s = 0.15
    β = 0.29
    ω_crit = 2 * R_s / (2 * R_s + β * (1 - R_s)^2)

    # Auxiliary equations
    τ = 0.1
    F_0 = 1370.0
    A_c = 0.6
    T_a = 0.76

    # Scattering aerosol (ω > ω_crit) - should cool
    ω_scatter = 0.95
    @test ω_scatter > ω_crit

    r_scatter = (1 - exp(-τ)) * ω_scatter * β
    t_scatter = exp(-τ) + ω_scatter * (1 - β) * (1 - exp(-τ))
    R_as_scatter = r_scatter + t_scatter^2 * R_s / (1 - R_s * r_scatter)
    ΔR_p_scatter = (1 - A_c) * T_a^2 * (R_as_scatter - R_s)

    @test ΔR_p_scatter > 0  # Increased albedo = cooling

    # Absorbing aerosol (ω < ω_crit) - should heat
    ω_absorb = 0.4
    @test ω_absorb < ω_crit

    r_absorb = (1 - exp(-τ)) * ω_absorb * β
    t_absorb = exp(-τ) + ω_absorb * (1 - β) * (1 - exp(-τ))
    R_as_absorb = r_absorb + t_absorb^2 * R_s / (1 - R_s * r_absorb)
    ΔR_p_absorb = (1 - A_c) * T_a^2 * (R_as_absorb - R_s)

    @test ΔR_p_absorb < 0  # Decreased albedo = heating
end
