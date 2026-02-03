@testsnippet TimescalesSetup begin
    using Test
    using ModelingToolkit
    using NonlinearSolve
    using Aerosol
end

@testitem "GasDiffusionTimescale structure" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = GasDiffusionTimescale()
    @test sys isa System
    @test length(equations(sys)) == 1
end

@testitem "GasDiffusionTimescale values" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = GasDiffusionTimescale()
    compiled = mtkcompile(sys)

    # Test τ_dg = R_p²/(4D_g)
    R_p = 1.0e-5  # 10 μm
    D_g = 2.0e-5  # m²/s

    expected_τ = R_p^2 / (4 * D_g)

    prob = NonlinearProblem(compiled, Dict(compiled.R_p => R_p, compiled.D_g => D_g))
    sol = solve(prob)

    @test sol[compiled.τ_dg] ≈ expected_τ rtol=1e-10
    # For 10 μm droplet, τ_dg ≈ 1.25e-6 s (very fast)
    @test sol[compiled.τ_dg] ≈ 1.25e-6 rtol=1e-6
end

@testitem "AqueousDiffusionTimescale structure" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = AqueousDiffusionTimescale()
    @test sys isa System
    @test length(equations(sys)) == 1
end

@testitem "AqueousDiffusionTimescale values" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = AqueousDiffusionTimescale()
    compiled = mtkcompile(sys)

    # Test τ_da = R_p²/(π²D_aq)
    R_p = 1.0e-5  # 10 μm
    D_aq = 1.0e-9  # m²/s (typical aqueous diffusivity)

    expected_τ = R_p^2 / (π^2 * D_aq)

    prob = NonlinearProblem(compiled, Dict(compiled.R_p => R_p, compiled.D_aq => D_aq))
    sol = solve(prob)

    @test sol[compiled.τ_da] ≈ expected_τ rtol=1e-10
    # For 10 μm droplet, τ_da ≈ 0.01 s
    @test sol[compiled.τ_da] ≈ 0.01 rtol=0.02
end

@testitem "InterfacialTimescale structure" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = InterfacialTimescale()
    @test sys isa System
    @test length(equations(sys)) == 2
end

@testitem "InterfacialTimescale ordering" setup=[TimescalesSetup] tags=[:timescales] begin
    # For most gases, τ_p_soluble > τ_p_insoluble when H* is large
    sys = InterfacialTimescale()
    compiled = mtkcompile(sys)

    # High Henry's law coefficient (very soluble)
    prob = NonlinearProblem(compiled, Dict(compiled.R_p => 1.0e-5, compiled.α => 1.0,
         compiled.M_A => 0.064, compiled.T => 298.15,
         compiled.H_star => 1.0e5, compiled.D_aq => 1.0e-9))
    sol = solve(prob)

    # Both timescales should be positive
    @test sol[compiled.τ_p_soluble] > 0
    @test sol[compiled.τ_p_insoluble] > 0
end

@testitem "ReactionTimescale structure" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = ReactionTimescale()
    @test sys isa System
    @test length(equations(sys)) == 2
end

@testitem "ReactionTimescale values" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = ReactionTimescale()
    compiled = mtkcompile(sys)

    # Test τ_ra = 1/k and τ_rg = τ_ra/(H*RT)
    k_rxn = 10.0  # s⁻¹
    H_star = 1.0e5  # mol/(m³·Pa)
    T = 298.15  # K
    R = 8.314  # J/(mol·K)

    expected_τ_ra = 1 / k_rxn
    expected_τ_rg = expected_τ_ra / (H_star * R * T)

    prob = NonlinearProblem(compiled, Dict(compiled.k_rxn => k_rxn, compiled.H_star => H_star, compiled.T => T))
    sol = solve(prob)

    @test sol[compiled.τ_ra] ≈ expected_τ_ra rtol=1e-10
    @test sol[compiled.τ_rg] ≈ expected_τ_rg rtol=1e-6
end

@testitem "SolidEquilibrationTimescale structure" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = SolidEquilibrationTimescale()
    @test sys isa System
    @test length(equations(sys)) == 2
end

@testitem "SolidEquilibrationTimescale values" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = SolidEquilibrationTimescale()
    compiled = mtkcompile(sys)

    # Test Eq. 12.135: τ_s = (ρ_p R_p²)/(3D_A m_p f(Kn,α))
    R_p = 1.0e-7  # 100 nm
    ρ_p = 1500.0  # kg/m³
    D_A = 2.0e-5  # m²/s
    m_p = 1.0e-8  # kg/m³
    f_Kn = 0.5  # transition regime correction

    expected_τ = (ρ_p * R_p^2) / (3 * D_A * m_p * f_Kn)

    prob = NonlinearProblem(compiled, Dict(compiled.R_p => R_p, compiled.ρ_p => ρ_p,
         compiled.D_A => D_A, compiled.m_p => m_p,
         compiled.f_Kn => f_Kn, compiled.N => 1.0e9))
    sol = solve(prob)

    @test sol[compiled.τ_s] ≈ expected_τ rtol=1e-10
end

@testitem "AqueousEquilibrationTimescale structure" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = AqueousEquilibrationTimescale()
    @test sys isa System
    @test length(equations(sys)) == 1
end

@testitem "AqueousEquilibrationTimescale values" setup=[TimescalesSetup] tags=[:timescales] begin
    sys = AqueousEquilibrationTimescale()
    compiled = mtkcompile(sys)

    # Test Eq. 12.147: τ_a = (m_w/K_A) τ_s
    m_w = 1.0e-8  # kg/m³
    K_A = 1.0e-6  # kg/m³
    τ_s = 1.0  # s

    expected_τ_a = (m_w / K_A) * τ_s

    prob = NonlinearProblem(compiled, Dict(compiled.m_w => m_w, compiled.K_A => K_A, compiled.τ_s => τ_s))
    sol = solve(prob)

    @test sol[compiled.τ_a] ≈ expected_τ_a rtol=1e-10
    # m_w/K_A = 1e-8/1e-6 = 0.01, so τ_a = 0.01 s
    @test sol[compiled.τ_a] ≈ 0.01 rtol=1e-6
end

@testitem "Timescale ordering" setup=[TimescalesSetup] tags=[:timescales] begin
    # For typical cloud droplet conditions, the ordering should be:
    # τ_dg < τ_da < τ_s (approximately)

    sys_dg = GasDiffusionTimescale()
    sys_da = AqueousDiffusionTimescale()

    compiled_dg = mtkcompile(sys_dg)
    compiled_da = mtkcompile(sys_da)

    R_p = 1.0e-5  # 10 μm droplet
    D_g = 2.0e-5
    D_aq = 1.0e-9

    prob_dg = NonlinearProblem(compiled_dg, Dict(compiled_dg.R_p => R_p, compiled_dg.D_g => D_g))
    sol_dg = solve(prob_dg)

    prob_da = NonlinearProblem(compiled_da, Dict(compiled_da.R_p => R_p, compiled_da.D_aq => D_aq))
    sol_da = solve(prob_da)

    # Gas-phase diffusion should be much faster than aqueous diffusion
    @test sol_dg[compiled_dg.τ_dg] < sol_da[compiled_da.τ_da]
    # The ratio should be approximately D_aq/D_g × (4/π²) ≈ 2e-5
    ratio = sol_dg[compiled_dg.τ_dg] / sol_da[compiled_da.τ_da]
    expected_ratio = (D_aq / D_g) * (π^2 / 4)
    @test ratio ≈ expected_ratio rtol=0.01
end
