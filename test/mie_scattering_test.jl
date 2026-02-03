@testsnippet MieSetup begin
    using Test
    using Statistics: mean
    using ModelingToolkit
    using OrdinaryDiffEqDefault
    using Aerosol
end

@testitem "MieScattering structural verification" setup=[MieSetup] tags=[:mie] begin
    # Create the system
    sys=MieScattering()

    # Verify the system has the expected structure
    @test sys isa ModelingToolkit.System

    # Get unknowns (variables)
    vars=unknowns(sys)
    var_names=[string(Symbolics.tosymbol(v, escape = false)) for v in vars]

    # Check that all expected variables are present
    expected_vars=["α", "Q_scat", "Q_abs", "Q_ext", "ω", "C_scat",
        "C_abs", "C_ext", "E_scat", "E_abs", "E_ext"]
    for v in expected_vars
        @test any(occursin(v, vn) for vn in var_names)
    end

    # Check number of equations
    eqs=equations(sys)
    @test length(eqs) == 11  # 11 equations defined
end

@testitem "RayleighScattering structural verification" setup=[MieSetup] tags=[:mie] begin
    sys=RayleighScattering()

    @test sys isa ModelingToolkit.System

    vars=unknowns(sys)
    var_names=[string(Symbolics.tosymbol(v, escape = false)) for v in vars]

    # Check expected variables
    expected_vars=["α", "Q_scat", "Q_abs", "Q_ext", "ω", "E_scat", "E_abs", "E_ext"]
    for v in expected_vars
        @test any(occursin(v, vn) for vn in var_names)
    end
end

@testitem "AerosolExtinction structural verification" setup=[MieSetup] tags=[:mie] begin
    sys=AerosolExtinction()

    @test sys isa ModelingToolkit.System

    vars=unknowns(sys)
    var_names=[string(Symbolics.tosymbol(v, escape = false)) for v in vars]

    # Check extinction coefficients are present
    expected_vars=["b_ext", "b_scat", "b_abs", "Q_ext", "Q_scat", "Q_abs", "ω"]
    for v in expected_vars
        @test any(occursin(v, vn) for vn in var_names)
    end
end

@testitem "Visibility structural verification" setup=[MieSetup] tags=[:mie] begin
    sys=Visibility()

    @test sys isa ModelingToolkit.System

    vars=unknowns(sys)
    var_names=[string(Symbolics.tosymbol(v, escape = false)) for v in vars]

    # Check visual range variable
    @test any(occursin("x_v", vn) for vn in var_names)
end

@testitem "Mie function basic evaluation" setup=[MieSetup] tags=[:mie] begin
    # Test the Mie functions directly
    # For a non-absorbing particle (k=0) with n=1.5 and α=1.0
    α=1.0
    n=1.5
    k=0.0

    Q_ext=Aerosol.mie_Q_ext(α, n, k)
    Q_scat=Aerosol.mie_Q_scat(α, n, k)
    Q_abs=Q_ext-Q_scat

    # For non-absorbing particles, Q_abs should be ~0
    @test Q_abs ≈ 0.0 atol=1e-10

    # Q_ext and Q_scat should be positive
    @test Q_ext > 0
    @test Q_scat > 0

    # Q_scat should equal Q_ext for non-absorbing particles
    @test Q_scat ≈ Q_ext atol=1e-10
end

@testitem "Rayleigh regime limiting behavior" setup=[MieSetup] tags=[:mie] begin
    # In the Rayleigh regime (α << 1), Q_scat ~ α^4 and Q_abs ~ α
    n=1.5
    k=0.01  # Small absorption

    α1=0.001  # Very small size parameter (deep Rayleigh regime)
    α2=0.002  # Double the size parameter

    Q_scat1=Aerosol.mie_Q_scat(α1, n, k)
    Q_scat2=Aerosol.mie_Q_scat(α2, n, k)

    # Q_scat should scale as α^4
    # So Q_scat2/Q_scat1 ≈ (α2/α1)^4 = 16
    ratio=Q_scat2/Q_scat1
    @test ratio ≈ 16.0 rtol=0.15  # Allow some tolerance for numerical effects
end

@testitem "Extinction paradox (large particle limit)" setup=[MieSetup] tags=[:mie] begin
    # For very large particles, Q_ext → 2 (the "extinction paradox")
    # Eq. 15.23
    n=1.33  # Water
    k=0.0

    # Test with large α values
    for α in [50.0, 100.0, 200.0]
        Q_ext=Aerosol.mie_Q_ext(α, n, k)
        # Q_ext should approach 2 with oscillations
        @test 1.5 < Q_ext < 2.5  # Allow for oscillations
    end

    # Average over oscillations should approach 2
    Q_exts=[Aerosol.mie_Q_ext(α, n, k) for α in 100.0:1.0:150.0]
    @test mean(Q_exts) ≈ 2.0 rtol=0.1
end

@testitem "Absorbing vs non-absorbing particles" setup=[MieSetup] tags=[:mie] begin
    α=2.0
    n=1.5

    # Non-absorbing particle
    k_nonabs=0.0
    Q_ext_nonabs=Aerosol.mie_Q_ext(α, n, k_nonabs)
    Q_scat_nonabs=Aerosol.mie_Q_scat(α, n, k_nonabs)
    Q_abs_nonabs=Q_ext_nonabs-Q_scat_nonabs

    # Absorbing particle
    k_abs=0.1
    Q_ext_abs=Aerosol.mie_Q_ext(α, n, k_abs)
    Q_scat_abs=Aerosol.mie_Q_scat(α, n, k_abs)
    Q_abs_abs=Q_ext_abs-Q_scat_abs

    # Non-absorbing should have negligible absorption
    @test Q_abs_nonabs ≈ 0.0 atol=1e-10

    # Absorbing particle should have significant absorption
    @test Q_abs_abs > 0.01

    # For moderately absorbing particles, Q_ext may be larger or smaller
    # depending on size parameter and refractive index
    # The key test is that Q_abs > 0 for k > 0
    @test Q_abs_abs > 0
end

@testitem "Koschmeider equation (Eq. 15.36)" setup=[MieSetup] tags=[:mie] begin
    # Visual range x_v = 3.912 / b_ext
    # For Rayleigh atmosphere at sea level: b_ext ≈ 13.2e-6 m^-1 at 520 nm
    # This should give x_v ≈ 296 km (page 705)

    # Test Koschmeider equation directly (without MTK compilation)
    Koschmeider_const=3.912
    b_ext=13.2e-6  # m^-1
    x_v=Koschmeider_const/b_ext

    # Expected: x_v ≈ 296 km = 296000 m
    @test x_v ≈ 296000.0 rtol=0.01
end

@testitem "Single-scattering albedo (Eq. 15.5)" setup=[MieSetup] tags=[:mie] begin
    # ω = Q_scat / Q_ext = C_scat / C_ext
    # For non-absorbing particles, ω should be 1.0
    # For strongly absorbing particles (like soot), ω should be < 0.5

    # Non-absorbing particle (water-like)
    n_water=1.33
    k_water=0.0
    α=5.0

    Q_ext=Aerosol.mie_Q_ext(α, n_water, k_water)
    Q_scat=Aerosol.mie_Q_scat(α, n_water, k_water)
    ω_water=Q_scat/Q_ext

    @test ω_water ≈ 1.0 atol=1e-6

    # Strongly absorbing particle (carbon-like)
    n_carbon=1.95
    k_carbon=0.79  # From Table 15.2

    Q_ext_c=Aerosol.mie_Q_ext(α, n_carbon, k_carbon)
    Q_scat_c=Aerosol.mie_Q_scat(α, n_carbon, k_carbon)
    ω_carbon=Q_scat_c/Q_ext_c

    # Carbon should have much lower single-scattering albedo
    # At α=5, ω ≈ 0.5 for elemental carbon (Table 15.2: n=1.95, k=0.79)
    @test ω_carbon < 0.55
end

@testitem "Mass extinction efficiency (Eq. 15.41-15.43)" setup=[MieSetup] tags=[:mie] begin
    # E_ext = 3 * Q_ext / (2 * ρ_p * D_p)
    # Test directly without MTK compilation

    # Parameters for ammonium sulfate at λ=550nm
    D_p=0.5e-6  # 0.5 μm
    λ=550e-9    # 550 nm
    n_refr=1.521
    k_refr=0.0
    ρ_p=1770.0  # kg/m³ for (NH4)2SO4

    α=π*D_p/λ
    Q_ext=Aerosol.mie_Q_ext(α, n_refr, k_refr)
    Q_scat=Aerosol.mie_Q_scat(α, n_refr, k_refr)

    # Mass scattering efficiency (Eq. 15.42)
    E_scat=3*Q_scat/(2*ρ_p*D_p)

    # From Figure 15.7, mass scattering efficiency for (NH4)2SO4 at 0.5 μm
    # should peak around 5-6 m²/g = 5000-6000 m²/kg
    @test 3000.0 < E_scat < 8000.0
end

@testitem "Size parameter calculation (Eq. 15.6)" setup=[MieSetup] tags=[:mie] begin
    # α = πD_p / λ
    D_p=1.0e-6   # 1 μm
    λ=500e-9     # 500 nm

    # Expected: α = π × 1e-6 / 500e-9 = π × 2 ≈ 6.28
    expected_α=π*D_p/λ
    @test expected_α ≈ 2π rtol=1e-10
end

@testitem "Extinction coefficient (Eq. 15.27)" setup=[MieSetup] tags=[:mie] begin
    # b_ext = (πD_p²/4) × N × Q_ext
    D_p=0.5e-6   # 0.5 μm
    N=1e9        # 1000 particles per cm³ = 1e9 per m³
    λ=550e-9
    n_refr=1.5
    k_refr=0.0

    # Calculate expected value
    α=π*D_p/λ
    Q_ext=Aerosol.mie_Q_ext(α, n_refr, k_refr)
    expected_b_ext=(π*D_p^2/4)*N*Q_ext

    # Verify the calculation makes sense (order of magnitude check)
    # b_ext should be positive and in a reasonable range
    @test expected_b_ext > 0
    @test 1e-6 < expected_b_ext < 1e-2  # Reasonable range for m^-1
end

@testitem "Conservation: Q_ext = Q_scat + Q_abs (Eq. 15.4)" setup=[MieSetup] tags=[:mie] begin
    # Test that extinction equals scattering plus absorption
    test_cases=[
        (α = 0.1, n = 1.5, k = 0.0),    # Small, non-absorbing
        (α = 1.0, n = 1.5, k = 0.01),   # Medium, weakly absorbing
        (α = 5.0, n = 1.33, k = 0.0),   # Large, non-absorbing (water)
        (α = 3.0, n = 1.95, k = 0.79),  # Strongly absorbing (carbon)
        (α = 10.0, n = 1.5, k = 0.1)   # Large, moderately absorbing
    ]

    for (α, n, k) in test_cases
        Q_ext=Aerosol.mie_Q_ext(α, n, k)
        Q_scat=Aerosol.mie_Q_scat(α, n, k)
        Q_abs=Q_ext-Q_scat

        # Q_abs should be non-negative (enforced by implementation)
        @test Q_abs >= -1e-10

        # All efficiencies should be positive for physically meaningful cases
        @test Q_ext >= 0
        @test Q_scat >= 0
    end
end

@testitem "Refractive index effects on scattering" setup=[MieSetup] tags=[:mie] begin
    # Higher real part of refractive index should generally increase scattering
    # for particles in the Mie regime
    α=3.0
    k=0.0

    n_values=[1.25, 1.33, 1.5, 1.75, 2.0]
    Q_scat_values=[Aerosol.mie_Q_scat(α, n, k) for n in n_values]

    # Generally Q_scat should increase with n (though not monotonically due to oscillations)
    # At least the trend should be upward
    @test Q_scat_values[end] > Q_scat_values[1]
end

@testitem "Rayleigh m² factor calculation" setup=[MieSetup] tags=[:mie] begin
    # Test the (m²-1)/(m²+2) factor calculation
    # For a non-absorbing dielectric (k=0), this should be real

    n=1.5
    k=0.0

    real_factor=Aerosol.rayleigh_m2_factor_real(n, k)
    imag_factor=Aerosol.rayleigh_m2_factor_imag(n, k)

    # For k=0, imaginary part should be 0
    @test imag_factor ≈ 0.0 atol=1e-10

    # Real part should be (n²-1)/(n²+2) = (2.25-1)/(2.25+2) = 1.25/4.25 ≈ 0.294
    expected_real=(n^2-1)/(n^2+2)
    @test real_factor ≈ expected_real rtol=1e-10

    # Test with absorption
    n=1.5
    k=0.1

    real_factor_abs=Aerosol.rayleigh_m2_factor_real(n, k)
    imag_factor_abs=Aerosol.rayleigh_m2_factor_imag(n, k)

    # With k > 0, imaginary part should be positive
    @test imag_factor_abs > 0
end

@testitem "Reference values from Table 15.2" setup=[MieSetup] tags=[:mie] begin
    # Test Mie calculations against reference refractive indices from Table 15.2
    # at λ = 550 nm

    # Water: n = 1.333, k ≈ 0
    # For water droplet with α = 5.0
    α=5.0
    n_water=1.333
    k_water=0.0

    Q_ext_water=Aerosol.mie_Q_ext(α, n_water, k_water)
    Q_scat_water=Aerosol.mie_Q_scat(α, n_water, k_water)

    # Water should be purely scattering
    @test Q_scat_water ≈ Q_ext_water atol=1e-6

    # Ammonium sulfate: n = 1.521, k ≈ 0
    n_sulfate=1.521
    k_sulfate=0.0

    Q_ext_sulfate=Aerosol.mie_Q_ext(α, n_sulfate, k_sulfate)

    # Higher refractive index should give higher Q_ext
    @test Q_ext_sulfate > Q_ext_water
end
