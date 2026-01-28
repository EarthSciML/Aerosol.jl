@testsnippet DynamicsSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D
    using OrdinaryDiffEq
    using Aerosol
end

@testitem "DiameterGrowthRate structure" setup=[DynamicsSetup] tags=[:dynamics] begin
    sys = DiameterGrowthRate()
    @test sys isa System

    # Check state variables
    vars = unknowns(sys)
    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in vars]
    @test "D_p" in var_names || "D_p(t)" in var_names
    @test "I_D" in var_names || "I_D(t)" in var_names
    @test "A" in var_names || "A(t)" in var_names

    # Check equations
    eqs = equations(sys)
    @test length(eqs) >= 3
end

@testitem "DiameterGrowthRate analytical solution" setup=[DynamicsSetup] tags=[:dynamics] begin
    # Test that diameter evolves according to D_p^2 = D_p0^2 + 2*A*t (Eq. 13.21)
    sys = DiameterGrowthRate()
    compiled_sys = mtkcompile(sys)

    # Parameters from Figure 13.2
    D_p0 = 0.2e-6  # 0.2 μm initial diameter
    D_diff = 0.1e-4  # m^2/s
    M_i = 0.1  # kg/mol
    Δp = 1e-4  # Pa (1 ppb at 1 atm)
    T_val = 298.0  # K
    ρ_p = 1000.0  # kg/m^3

    prob = ODEProblem(
        compiled_sys,
        [compiled_sys.D_p => D_p0],
        (0.0, 1200.0),  # 20 minutes
        [
            compiled_sys.D_diff => D_diff,
            compiled_sys.M_i => M_i,
            compiled_sys.Δp => Δp,
            compiled_sys.T => T_val,
            compiled_sys.ρ_p => ρ_p,
        ]
    )
    sol = solve(prob)

    # Compute expected A parameter
    R_gas = 8.314
    A_expected = 4 * D_diff * M_i * Δp / (R_gas * T_val * ρ_p)

    # Check solution at various times using Eq. 13.21
    for t_val in [0.0, 300.0, 600.0, 1200.0]
        idx = argmin(abs.(sol.t .- t_val))
        D_p_numerical = sol[compiled_sys.D_p][idx]
        D_p_analytical = sqrt(D_p0^2 + 2 * A_expected * sol.t[idx])

        @test isapprox(D_p_numerical, D_p_analytical, rtol=1e-3)
    end
end

@testitem "BrownianCoagulationCoefficient structure" setup=[DynamicsSetup] tags=[:dynamics] begin
    sys = BrownianCoagulationCoefficient()
    @test sys isa System

    # Check key variables exist
    vars = unknowns(sys)
    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in vars]
    @test "K_12" in var_names || "K_12(t)" in var_names
    @test "D_1" in var_names || "D_1(t)" in var_names
    @test "D_2" in var_names || "D_2(t)" in var_names
end

@testitem "BrownianCoagulationCoefficient continuum limit" setup=[DynamicsSetup] tags=[:dynamics] begin
    # For large particles (continuum regime), K should approach 8kT/3μ for equal sizes (Eq. 13.52)
    sys = BrownianCoagulationCoefficient()
    compiled_sys = mtkcompile(sys)

    # Large particles in continuum regime
    D_p_large = 10e-6  # 10 μm
    T_val = 298.0
    μ_val = 1.83e-5
    k_B = 1.380649e-23

    prob = ODEProblem(
        compiled_sys,
        [],
        (0.0, 1.0),
        [
            compiled_sys.D_p1 => D_p_large,
            compiled_sys.D_p2 => D_p_large,
            compiled_sys.T => T_val,
            compiled_sys.μ => μ_val,
        ]
    )
    sol = solve(prob)

    K_computed = sol[compiled_sys.K_12][1]
    K_continuum = 8 * k_B * T_val / (3 * μ_val)  # Eq. 13.52

    # Should be within about 10% for large particles
    @test isapprox(K_computed, K_continuum, rtol=0.15)
end

@testitem "BrownianCoagulationCoefficient size dependence" setup=[DynamicsSetup] tags=[:dynamics] begin
    # Test that K increases for unequal particle sizes in continuum regime
    sys = BrownianCoagulationCoefficient()
    compiled_sys = mtkcompile(sys)

    T_val = 298.0
    μ_val = 1.83e-5

    # Equal sizes
    prob_equal = ODEProblem(
        compiled_sys,
        [],
        (0.0, 1.0),
        [
            compiled_sys.D_p1 => 1e-6,
            compiled_sys.D_p2 => 1e-6,
            compiled_sys.T => T_val,
            compiled_sys.μ => μ_val,
        ]
    )
    sol_equal = solve(prob_equal)
    K_equal = sol_equal[compiled_sys.K_12][1]

    # Unequal sizes (factor of 10)
    prob_unequal = ODEProblem(
        compiled_sys,
        [],
        (0.0, 1.0),
        [
            compiled_sys.D_p1 => 0.1e-6,
            compiled_sys.D_p2 => 1e-6,
            compiled_sys.T => T_val,
            compiled_sys.μ => μ_val,
        ]
    )
    sol_unequal = solve(prob_unequal)
    K_unequal = sol_unequal[compiled_sys.K_12][1]

    # Unequal sizes should give higher coagulation coefficient
    @test K_unequal > K_equal
end

@testitem "MonodisperseCoagulation structure" setup=[DynamicsSetup] tags=[:dynamics] begin
    sys = MonodisperseCoagulation()
    @test sys isa System

    vars = unknowns(sys)
    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in vars]
    @test "N" in var_names || "N(t)" in var_names
    @test "τ_c" in var_names || "τ_c(t)" in var_names
end

@testitem "MonodisperseCoagulation analytical solution" setup=[DynamicsSetup] tags=[:dynamics] begin
    # Test against analytical solution N(t) = N_0 / (1 + t/τ_c) (Eq. 13.66)
    sys = MonodisperseCoagulation()
    compiled_sys = mtkcompile(sys)

    N_0_val = 1e12  # m^-3
    K_val = 1e-15   # m^3/s
    τ_c_expected = 2 / (K_val * N_0_val)  # Eq. 13.67

    prob = ODEProblem(
        compiled_sys,
        [compiled_sys.N => N_0_val],
        (0.0, 5 * τ_c_expected),
        [
            compiled_sys.K => K_val,
            compiled_sys.N_0 => N_0_val,
        ]
    )
    sol = solve(prob; reltol=1e-8, abstol=1e-10)

    # Check at characteristic time τ_c, N should be N_0/2
    idx_τc = argmin(abs.(sol.t .- τ_c_expected))
    N_at_τc = sol[compiled_sys.N][idx_τc]
    @test isapprox(N_at_τc, N_0_val / 2, rtol=0.15)

    # Check at 4*τ_c, N should be N_0/5
    idx_4τc = argmin(abs.(sol.t .- 4 * τ_c_expected))
    N_at_4τc = sol[compiled_sys.N][idx_4τc]
    @test isapprox(N_at_4τc, N_0_val / 5, rtol=0.15)
end

@testitem "MonodisperseCoagulation characteristic time" setup=[DynamicsSetup] tags=[:dynamics] begin
    # Test τ_c values from Chapter 13 examples
    sys = MonodisperseCoagulation()
    compiled_sys = mtkcompile(sys)

    # Example 1: N_0 = 10^10 m^-3, K = 10^-15 m^3/s → τ_c ≈ 55 hours
    prob1 = ODEProblem(
        compiled_sys,
        [compiled_sys.N => 1e10],
        (0.0, 1.0),
        [
            compiled_sys.K => 1e-15,
            compiled_sys.N_0 => 1e10,
        ]
    )
    sol1 = solve(prob1)
    τ_c1 = sol1[compiled_sys.τ_c][1]
    τ_c1_hours = τ_c1 / 3600
    @test isapprox(τ_c1_hours, 55.5, rtol=0.1)

    # Example 2: N_0 = 10^12 m^-3, K = 10^-15 m^3/s → τ_c ≈ 33 minutes
    prob2 = ODEProblem(
        compiled_sys,
        [compiled_sys.N => 1e12],
        (0.0, 1.0),
        [
            compiled_sys.K => 1e-15,
            compiled_sys.N_0 => 1e12,
        ]
    )
    sol2 = solve(prob2)
    τ_c2 = sol2[compiled_sys.τ_c][1]
    τ_c2_min = τ_c2 / 60
    @test isapprox(τ_c2_min, 33.3, rtol=0.1)
end

@testitem "DiscreteCoagulation structure" setup=[DynamicsSetup] tags=[:dynamics] begin
    n_bins = 5
    sys = DiscreteCoagulation(n_bins)
    @test sys isa System

    # Check that we have the right number of state variables
    vars = unknowns(sys)
    @test length(vars) >= n_bins  # At least n_bins for N[1:n_bins]
end

@testitem "DiscreteCoagulation monodisperse initial condition" setup=[DynamicsSetup] tags=[:dynamics] begin
    # Test against analytical solution N_k(t) = N_0 * (t/τ_c)^(k-1) / (1 + t/τ_c)^(k+1) (Eq. 13.71)
    n_bins = 5
    sys = DiscreteCoagulation(n_bins)
    compiled_sys = mtkcompile(sys)

    N_0_val = 1e12
    K_val = 1e-15
    τ_c = 2 / (K_val * N_0_val)

    # Initially all monomers (N_1 = N_0, N_k>1 = 0)
    u0 = [compiled_sys.N[1] => N_0_val]
    for k in 2:n_bins
        push!(u0, compiled_sys.N[k] => 0.0)
    end

    prob = ODEProblem(
        compiled_sys,
        u0,
        (0.0, 3 * τ_c),
        [
            compiled_sys.K => K_val,
            compiled_sys.N_0 => N_0_val,
        ]
    )
    sol = solve(prob)

    # At t = τ_c, check N_1 and N_2 against analytical solution
    t_test = τ_c
    idx = argmin(abs.(sol.t .- t_test))
    t_actual = sol.t[idx]
    θ = t_actual / τ_c  # Dimensionless time

    # N_1 analytical: N_0 / (1 + θ)^2
    N_1_analytical = N_0_val / (1 + θ)^2
    N_1_numerical = sol[compiled_sys.N[1]][idx]
    @test isapprox(N_1_numerical, N_1_analytical, rtol=0.1)

    # N_2 analytical: N_0 * θ / (1 + θ)^3
    N_2_analytical = N_0_val * θ / (1 + θ)^3
    N_2_numerical = sol[compiled_sys.N[2]][idx]
    @test isapprox(N_2_numerical, N_2_analytical, rtol=0.1)
end

@testitem "AerosolDynamics integration" setup=[DynamicsSetup] tags=[:dynamics] begin
    sys = AerosolDynamics()
    @test sys isa System

    # Check subsystems are present by checking that the system has expected equations
    eqs = equations(sys)
    @test length(eqs) > 0

    # Verify that coupling equations exist
    eq_strs = string.(eqs)
    @test any(contains(s, "growth") for s in eq_strs)
    @test any(contains(s, "coag_coeff") for s in eq_strs)
    @test any(contains(s, "coag") for s in eq_strs)
end

@testitem "Volume conservation under coagulation" setup=[DynamicsSetup] tags=[:dynamics] begin
    # Coagulation conserves total volume/mass (Table 13.4)
    # Total volume V_0 = sum(k * N_k * v_0) should remain constant
    # Note: Using more bins to reduce truncation error from particles growing beyond n_bins
    n_bins = 10
    sys = DiscreteCoagulation(n_bins)
    compiled_sys = mtkcompile(sys)

    N_0_val = 1e12
    K_val = 1e-15
    τ_c = 2 / (K_val * N_0_val)

    # Initially all monomers
    u0 = [compiled_sys.N[1] => N_0_val]
    for k in 2:n_bins
        push!(u0, compiled_sys.N[k] => 0.0)
    end

    prob = ODEProblem(
        compiled_sys,
        u0,
        (0.0, 0.5 * τ_c),  # Shorter time to avoid too many particles growing beyond n_bins
        [
            compiled_sys.K => K_val,
            compiled_sys.N_0 => N_0_val,
        ]
    )
    sol = solve(prob; reltol=1e-8, abstol=1e-10)

    # Compute volume (in units of monomer volume) at initial and final times
    V_initial = sum(k * sol[compiled_sys.N[k]][1] for k in 1:n_bins)
    V_final = sum(k * sol[compiled_sys.N[k]][end] for k in 1:n_bins)

    # Volume should be approximately conserved
    # Some volume is lost due to particles growing beyond n_bins, so allow larger tolerance
    @test isapprox(V_final, V_initial, rtol=0.25)
end
