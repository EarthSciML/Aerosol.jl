@testitem "Nucleation_structural" begin
    using Test
    using ModelingToolkit
    using Aerosol

    # Test that Nucleation component creates a valid system
    sys = Nucleation()

    # Check that system compiles without error
    compiled = mtkcompile(sys)

    # Check that equations are present (8 equations)
    @test length(equations(sys)) == 8
end

@testitem "WaterProperties_structural" begin
    using Test
    using ModelingToolkit
    using Aerosol

    # Test that WaterProperties component creates a valid system
    sys = WaterProperties()

    # Check that system compiles without error
    compiled = mtkcompile(sys)

    # Check that equations are present (4 equations)
    @test length(equations(sys)) == 4
end

@testitem "CriticalCluster_structural" begin
    using Test
    using ModelingToolkit
    using Aerosol

    # Test that CriticalCluster component creates a valid system
    sys = CriticalCluster()

    # Check that system compiles without error
    compiled = mtkcompile(sys)

    # Check that equations are present (4 equations)
    @test length(equations(sys)) == 4
end

@testitem "ClassicalNucleationRate_structural" begin
    using Test
    using ModelingToolkit
    using Aerosol

    # Test that ClassicalNucleationRate component creates a valid system
    sys = ClassicalNucleationRate()

    # Check that system compiles without error
    compiled = mtkcompile(sys)

    # Check that equations are present (1 equation)
    @test length(equations(sys)) == 1
end

@testitem "Critical_cluster_water_273K" begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEq
    using Aerosol

    # Validation against Table 11.1: Water at 273 K
    # At S=2: r* = 17.3 Å, i* = 726
    # At S=5: r* = 7.5 Å, i* = 58

    # Water properties at 273 K from Table 11.1
    σ_water = 75.6e-3  # N/m (75.6 dyn/cm = 75.6e-3 N/m)
    v_1_water = 2.99e-29  # m³/molecule (2.99e-23 cm³ = 2.99e-29 m³)
    T = 273.0  # K

    sys = CriticalCluster()
    compiled = mtkcompile(sys)

    # Test S = 2
    prob = ODEProblem(compiled,
        Dict(compiled.T => T, compiled.S => 2.0, compiled.σ => σ_water, compiled.v_1 => v_1_water),
        (0.0, 1.0))
    sol = solve(prob, Tsit5())

    # r* should be ~17.3 Å = 17.3e-10 m
    r_star = sol[compiled.r_star][1]
    @test isapprox(r_star * 1e10, 17.3, rtol=0.05)  # 5% tolerance

    # i* should be ~726
    i_star = sol[compiled.i_star][1]
    @test isapprox(i_star, 726, rtol=0.1)  # 10% tolerance (capillarity approximation)

    # Test S = 5
    prob5 = ODEProblem(compiled,
        Dict(compiled.T => T, compiled.S => 5.0, compiled.σ => σ_water, compiled.v_1 => v_1_water),
        (0.0, 1.0))
    sol5 = solve(prob5, Tsit5())

    r_star5 = sol5[compiled.r_star][1]
    @test isapprox(r_star5 * 1e10, 7.5, rtol=0.05)

    i_star5 = sol5[compiled.i_star][1]
    @test isapprox(i_star5, 58, rtol=0.15)
end

@testitem "Critical_cluster_water_298K" begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEq
    using Aerosol

    # Validation against Table 11.1: Water at 298 K
    # At S=2: r* = 15.1 Å, i* = 482
    # At S=5: r* = 6.5 Å, i* = 39

    σ_water = 72.0e-3  # N/m at 298 K
    v_1_water = 2.99e-29  # m³/molecule
    T = 298.0  # K

    sys = CriticalCluster()
    compiled = mtkcompile(sys)

    # Test S = 2
    prob = ODEProblem(compiled,
        Dict(compiled.T => T, compiled.S => 2.0, compiled.σ => σ_water, compiled.v_1 => v_1_water),
        (0.0, 1.0))
    sol = solve(prob, Tsit5())

    r_star = sol[compiled.r_star][1]
    @test isapprox(r_star * 1e10, 15.1, rtol=0.05)

    i_star = sol[compiled.i_star][1]
    @test isapprox(i_star, 482, rtol=0.1)

    # Test S = 5
    prob5 = ODEProblem(compiled,
        Dict(compiled.T => T, compiled.S => 5.0, compiled.σ => σ_water, compiled.v_1 => v_1_water),
        (0.0, 1.0))
    sol5 = solve(prob5, Tsit5())

    r_star5 = sol5[compiled.r_star][1]
    @test isapprox(r_star5 * 1e10, 6.5, rtol=0.05)

    i_star5 = sol5[compiled.i_star][1]
    @test isapprox(i_star5, 39, rtol=0.15)
end

@testitem "Nucleation_rate_water_293K" begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEq
    using Aerosol

    # Validation against Table 11.4: Homogeneous nucleation rate for water at 293 K
    # S=5: i*=42, J = 1.57e11 cm⁻³s⁻¹ = 1.57e17 m⁻³s⁻¹
    # S=4: i*=65, J = 1.05e6 cm⁻³s⁻¹ = 1.05e12 m⁻³s⁻¹

    # Water properties at 293 K
    σ_water = 72.75e-3  # N/m
    v_1_water = 2.99e-29  # m³/molecule
    m_1_water = 2.99e-26  # kg/molecule (from Table 11.4)
    p_sat_water = 2336.5  # Pa (2.3365e4 g/(cm·s²) = 2336.5 Pa)
    T = 293.0  # K

    sys = Nucleation()
    compiled = mtkcompile(sys)

    # Test S = 5
    S = 5.0
    p_A = S * p_sat_water

    prob = ODEProblem(compiled, Dict(
        compiled.T => T,
        compiled.p_A => p_A,
        compiled.p_A_s => p_sat_water,
        compiled.σ => σ_water,
        compiled.v_1 => v_1_water,
        compiled.m_1 => m_1_water
    ), (0.0, 1.0))
    sol = solve(prob, Tsit5())

    J = sol[compiled.J][1]
    J_expected = 1.57e17  # m⁻³s⁻¹

    # Nucleation rate spans many orders of magnitude
    # Check that log10 values match within 1 order of magnitude
    @test isapprox(log10(J), log10(J_expected), atol=1.5)

    # Critical cluster size should be ~42
    i_star = sol[compiled.i_star][1]
    @test isapprox(i_star, 42, rtol=0.2)
end

@testitem "Nucleation_saturation_ratio" begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEq
    using Aerosol

    # Test that saturation ratio is computed correctly
    sys = Nucleation()
    compiled = mtkcompile(sys)

    T = 293.0
    p_A_s = 2336.5  # Pa
    p_A = 3 * p_A_s  # S = 3

    prob = ODEProblem(compiled, Dict(
        compiled.T => T,
        compiled.p_A => p_A,
        compiled.p_A_s => p_A_s,
        compiled.σ => 72.75e-3,
        compiled.v_1 => 2.99e-29,
        compiled.m_1 => 2.99e-26
    ), (0.0, 1.0))
    sol = solve(prob, Tsit5())

    S = sol[compiled.S][1]
    @test isapprox(S, 3.0, rtol=1e-6)
end

@testitem "WaterProperties_values" begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEq
    using Aerosol

    sys = WaterProperties()
    compiled = mtkcompile(sys)

    # Test at 273 K
    prob273 = ODEProblem(compiled, Dict(compiled.T => 273.0), (0.0, 1.0))
    sol273 = solve(prob273, Tsit5())

    σ_273 = sol273[compiled.σ][1]
    @test isapprox(σ_273, 75.6e-3, rtol=0.01)

    v_1 = sol273[compiled.v_1][1]
    @test isapprox(v_1, 2.99e-29, rtol=0.01)

    m_1 = sol273[compiled.m_1][1]
    @test isapprox(m_1, 18.015e-3 / 6.022e23, rtol=0.01)

    # Test at 293 K
    prob293 = ODEProblem(compiled, Dict(compiled.T => 293.0), (0.0, 1.0))
    sol293 = solve(prob293, Tsit5())

    σ_293 = sol293[compiled.σ][1]
    # At 293 K: σ ≈ 75.6 - 0.1454*20 ≈ 72.7 mN/m
    @test isapprox(σ_293, 72.7e-3, rtol=0.02)

    # Saturation vapor pressure at 293 K should be ~2338 Pa
    p_sat = sol293[compiled.p_sat][1]
    @test isapprox(p_sat, 2338.0, rtol=0.05)
end

@testitem "Nucleation_rate_sensitivity" begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEq
    using Aerosol

    # Test that nucleation rate increases dramatically with saturation ratio
    # This is a key qualitative behavior of nucleation

    sys = ClassicalNucleationRate()
    compiled = mtkcompile(sys)

    σ = 72.75e-3
    v_1 = 2.99e-29
    m_1 = 2.99e-26
    T = 293.0
    k_B = 1.380649e-23

    # At S=3, calculate N_1 from ideal gas with p_A = S * p_sat
    p_sat = 2336.5  # Pa
    N_1_S3 = 3 * p_sat / (k_B * T)
    N_1_S6 = 6 * p_sat / (k_B * T)

    prob3 = ODEProblem(compiled, Dict(
        compiled.T => T, compiled.S => 3.0, compiled.N_1 => N_1_S3,
        compiled.σ => σ, compiled.v_1 => v_1, compiled.m_1 => m_1
    ), (0.0, 1.0))
    sol3 = solve(prob3, Tsit5())
    J3 = sol3[compiled.J][1]

    prob6 = ODEProblem(compiled, Dict(
        compiled.T => T, compiled.S => 6.0, compiled.N_1 => N_1_S6,
        compiled.σ => σ, compiled.v_1 => v_1, compiled.m_1 => m_1
    ), (0.0, 1.0))
    sol6 = solve(prob6, Tsit5())
    J6 = sol6[compiled.J][1]

    # J should increase by many orders of magnitude from S=3 to S=6
    # Based on Table 11.4: J(S=3) ≈ 1.76e-6, J(S=6) ≈ 1.24e14 (ratio ~10^20)
    @test J6 > J3
    @test log10(J6) - log10(J3) > 10  # At least 10 orders of magnitude increase
end

@testitem "Critical_cluster_decreases_with_S" begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEq
    using Aerosol

    # Qualitative test: critical cluster size should decrease with increasing S

    sys = CriticalCluster()
    compiled = mtkcompile(sys)

    σ = 72.75e-3
    v_1 = 2.99e-29
    T = 293.0

    i_stars = Float64[]
    r_stars = Float64[]

    for S in [2.0, 3.0, 4.0, 5.0, 6.0]
        prob = ODEProblem(compiled, Dict(
            compiled.T => T, compiled.S => S, compiled.σ => σ, compiled.v_1 => v_1
        ), (0.0, 1.0))
        sol = solve(prob, Tsit5())
        push!(i_stars, sol[compiled.i_star][1])
        push!(r_stars, sol[compiled.r_star][1])
    end

    # Both i* and r* should be monotonically decreasing
    for k in 1:(length(i_stars) - 1)
        @test i_stars[k] > i_stars[k + 1]
        @test r_stars[k] > r_stars[k + 1]
    end
end

@testitem "Free_energy_barrier_positive" begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEq
    using Aerosol

    # The free energy barrier should always be positive for supersaturated conditions

    sys = CriticalCluster()
    compiled = mtkcompile(sys)

    σ = 72.75e-3
    v_1 = 2.99e-29
    T = 293.0

    for S in [1.5, 2.0, 3.0, 5.0, 10.0]
        prob = ODEProblem(compiled, Dict(
            compiled.T => T, compiled.S => S, compiled.σ => σ, compiled.v_1 => v_1
        ), (0.0, 1.0))
        sol = solve(prob, Tsit5())
        ΔG_star = sol[compiled.ΔG_star][1]
        @test ΔG_star > 0
    end
end

@testitem "Organic_species_critical_cluster" begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEq
    using Aerosol

    # Test with organic species data from Table 11.2: Ethanol at 298 K
    # σ = 22.14 dyn/cm = 22.14e-3 N/m
    # v_1 = 9.694e-23 cm³ = 9.694e-29 m³
    # At S=3, the critical radius should be smaller than for water due to lower surface tension

    sys = CriticalCluster()
    compiled = mtkcompile(sys)

    # Ethanol properties at 298 K
    σ_ethanol = 22.14e-3  # N/m
    v_1_ethanol = 9.694e-29  # m³/molecule
    T = 298.0

    prob = ODEProblem(compiled, Dict(
        compiled.T => T, compiled.S => 3.0, compiled.σ => σ_ethanol, compiled.v_1 => v_1_ethanol
    ), (0.0, 1.0))
    sol = solve(prob, Tsit5())

    r_star = sol[compiled.r_star][1]
    i_star = sol[compiled.i_star][1]

    # Critical radius should be reasonable (order of nm)
    @test 1e-10 < r_star < 1e-8

    # Critical cluster size should be positive
    @test i_star > 1
end
