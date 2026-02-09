@testsnippet DahnekeSetup begin
    using Test
    using ModelingToolkit
    using NonlinearSolve
    using Aerosol
end

# ============================================================
# DahnekeMassTransportCorrection tests
# ============================================================

@testitem "DahnekeMassTransportCorrection structure" setup=[DahnekeSetup] tags=[:dahneke] begin
    sys = DahnekeMassTransportCorrection()
    @test sys isa System
    @test length(equations(sys)) == 4  # c_bar, ℓ_D, Kn_D, β
end

@testitem "DahnekeMassTransportCorrection continuum limit" setup=[DahnekeSetup] tags=[:dahneke] begin
    # Eq. 5.5: β → 1 as Kn_D → 0
    sys = DahnekeMassTransportCorrection()
    compiled = mtkcompile(sys)

    # Use large sphere radius so Kn_D ≈ 0
    # For D_v=2e-5, T=298 K, m_v=4.81e-26 kg:
    # c_bar ≈ 417 m/s, ℓ_D = 2D/c̄ ≈ 9.6e-8 m
    # For r = 0.01 m: Kn_D = ℓ_D/r ≈ 1e-5
    prob = NonlinearProblem(compiled,
        Dict(compiled.r => 0.01, compiled.δ_m => 1.0))
    sol = solve(prob)
    @test sol[compiled.β] ≈ 1.0 rtol = 1e-4
    @test sol[compiled.Kn_D] < 1e-4
end

@testitem "DahnekeMassTransportCorrection free-molecular limit" setup=[DahnekeSetup] tags=[:dahneke] begin
    # Eq. 5.5: as Kn_D → ∞, β → δ/(2Kn_D) → δc̄r/(4D)
    sys = DahnekeMassTransportCorrection()
    compiled = mtkcompile(sys)

    # Use very small sphere so Kn_D >> 1
    prob = NonlinearProblem(compiled,
        Dict(compiled.r => 1e-9, compiled.δ_m => 1.0))
    sol = solve(prob)
    Kn_D_val = sol[compiled.Kn_D]
    β_val = sol[compiled.β]
    @test Kn_D_val > 10  # Should be in free-molecular regime

    # In free-molecular limit: β ≈ (Kn_D+1)/(2Kn_D²/δ + 1) ≈ δ/(2Kn_D)
    β_expected = 1.0 / (2 * Kn_D_val)
    @test β_val ≈ β_expected rtol = 0.05
end

@testitem "DahnekeMassTransportCorrection δ dependence" setup=[DahnekeSetup] tags=[:dahneke] begin
    # Lower sticking coefficient should reduce β
    sys = DahnekeMassTransportCorrection()
    compiled = mtkcompile(sys)

    # Choose r so Kn_D is moderate (transition regime)
    prob_δ1 = NonlinearProblem(compiled,
        Dict(compiled.r => 1e-7, compiled.δ_m => 1.0))
    sol_δ1 = solve(prob_δ1)

    prob_δ01 = NonlinearProblem(compiled,
        Dict(compiled.r => 1e-7, compiled.δ_m => 0.1))
    sol_δ01 = solve(prob_δ01)

    @test sol_δ01[compiled.β] < sol_δ1[compiled.β]
end

@testitem "DahnekeMassTransportCorrection exact formula" setup=[DahnekeSetup] tags=[:dahneke] begin
    # Verify Eq. 5.5 exactly: β = (Kn_D + 1) / (2*Kn_D*(Kn_D+1)/δ + 1)
    sys = DahnekeMassTransportCorrection()
    compiled = mtkcompile(sys)

    for (r_val, δ_val) in [(1e-7, 1.0), (1e-8, 0.5), (1e-6, 0.01)]
        prob = NonlinearProblem(compiled,
            Dict(compiled.r => r_val, compiled.δ_m => δ_val))
        sol = solve(prob)
        Kn = sol[compiled.Kn_D]
        β_calc = (Kn + 1) / (2 * Kn * (Kn + 1) / δ_val + 1)
        @test sol[compiled.β] ≈ β_calc rtol = 1e-10
    end
end

# ============================================================
# DahnekeHeatTransportCorrection tests
# ============================================================

@testitem "DahnekeHeatTransportCorrection structure" setup=[DahnekeSetup] tags=[:dahneke] begin
    sys = DahnekeHeatTransportCorrection()
    @test sys isa System
    @test length(equations(sys)) == 2  # Kn_K, β_q
end

@testitem "DahnekeHeatTransportCorrection exact formula" setup=[DahnekeSetup] tags=[:dahneke] begin
    # Verify Eq. 5.7: β_q = (Kn_K + 1) / (2*Kn_K*(Kn_K+1)/α + 1)
    sys = DahnekeHeatTransportCorrection()
    compiled = mtkcompile(sys)

    prob = NonlinearProblem(compiled,
        Dict(compiled.r => 1e-7, compiled.α => 1.0))
    sol = solve(prob)
    Kn_K = sol[compiled.Kn_K]
    β_q_calc = (Kn_K + 1) / (2 * Kn_K * (Kn_K + 1) / 1.0 + 1)
    @test sol[compiled.β_q] ≈ β_q_calc rtol = 1e-10
end

# ============================================================
# DahnekeCondensationEvaporation tests
# ============================================================

@testitem "DahnekeCondensationEvaporation structure" setup=[DahnekeSetup] tags=[:dahneke] begin
    sys = DahnekeCondensationEvaporation()
    @test sys isa System
end

@testitem "DahnekeCondensationEvaporation Maxwell equation" setup=[DahnekeSetup] tags=[:dahneke] begin
    # Verify I_M = 4πrD(n_∞ - n_s) (Eq. 2.1)
    sys = DahnekeCondensationEvaporation()
    compiled = mtkcompile(sys)

    r_val = 1e-6
    D_val = 2e-5
    n_inf_val = 1e21
    n_s_val = 0.0
    prob = NonlinearProblem(compiled,
        Dict(compiled.r => r_val, compiled.D_v => D_val,
            compiled.n_inf => n_inf_val, compiled.n_s => n_s_val))
    sol = solve(prob)

    I_M_expected = 4 * π * r_val * D_val * (n_inf_val - n_s_val)
    @test sol[compiled.I_M] ≈ I_M_expected rtol = 1e-8

    # Corrected rate should be I = β * I_M
    β = sol[compiled.mass_corr.β]
    @test sol[compiled.I] ≈ β * I_M_expected rtol = 1e-8
end

@testitem "DahnekeCondensationEvaporation heat equation" setup=[DahnekeSetup] tags=[:dahneke] begin
    # Verify Q_M = 4πrκ(T_∞ - T_o) (Eq. 2.2)
    sys = DahnekeCondensationEvaporation()
    compiled = mtkcompile(sys)

    r_val = 1e-6
    κ_val = 0.026
    T_inf_val = 298.15
    T_o_val = 293.15
    prob = NonlinearProblem(compiled,
        Dict(compiled.r => r_val, compiled.κ => κ_val,
            compiled.T_inf => T_inf_val, compiled.T_o => T_o_val))
    sol = solve(prob)

    Q_M_expected = 4 * π * r_val * κ_val * (T_inf_val - T_o_val)
    @test sol[compiled.Q_M] ≈ Q_M_expected rtol = 1e-8
end

# ============================================================
# DahnekeCoagulationRate tests
# ============================================================

@testitem "DahnekeCoagulationRate structure" setup=[DahnekeSetup] tags=[:dahneke] begin
    sys = DahnekeCoagulationRate()
    @test sys isa System

    # Check equation count
    eqs = equations(sys)
    @test length(eqs) == 16
end

@testitem "DahnekeCoagulationRate continuum limit" setup=[DahnekeSetup] tags=[:dahneke] begin
    # When Kn_D → 0, K → K_o = 4πRD (Eq. 8.8)
    sys = DahnekeCoagulationRate()
    compiled = mtkcompile(sys)

    # Large particles → small Kn_D
    prob = NonlinearProblem(compiled,
        Dict(compiled.r_1 => 1e-4, compiled.r_2 => 1e-4,
            compiled.C_s1 => 1.0, compiled.C_s2 => 1.0,
            compiled.δ_p => 1.0))
    sol = solve(prob)

    @test sol[compiled.Kn_D] < 0.01
    @test sol[compiled.K] ≈ sol[compiled.K_o] rtol = 0.02
end

@testitem "DahnekeCoagulationRate equal-size formula" setup=[DahnekeSetup] tags=[:dahneke] begin
    # For equal-size particles with C_s1 = C_s2 = 1:
    # K_o = 8πrD = 8kT/(3η) (Smoluchowski result when C_s=1)
    sys = DahnekeCoagulationRate()
    compiled = mtkcompile(sys)

    T_val = 298.15
    μ_val = 1.84e-5
    r_val = 1e-4  # Large particle for continuum limit

    prob = NonlinearProblem(compiled,
        Dict(compiled.T => T_val, compiled.μ => μ_val,
            compiled.r_1 => r_val, compiled.r_2 => r_val,
            compiled.C_s1 => 1.0, compiled.C_s2 => 1.0,
            compiled.δ_p => 1.0))
    sol = solve(prob)

    # K_o should equal Smoluchowski: K_o = 4πRD where R=2r and D=2*kT/(6πηr)
    k_B = 1.380649e-23
    R = 2 * r_val
    D_each = k_B * T_val / (6 * π * μ_val * r_val)
    D_mut = 2 * D_each
    K_o_expected = 4 * π * R * D_mut  # = 8kT/(3η)
    K_o_smoluchowski = 8 * k_B * T_val / (3 * μ_val)

    @test sol[compiled.K_o] ≈ K_o_expected rtol = 1e-6
    @test sol[compiled.K_o] ≈ K_o_smoluchowski rtol = 1e-6
end

@testitem "DahnekeCoagulationRate Table 2 validation" setup=[DahnekeSetup] tags=[:dahneke] begin
    # Validate against Table 2 from Dahneke 1983
    # Parameters: r₁=r₂=0.22 μm, T=25°C, δ=1, η=1.84e-5 Pa·s, ρ=917 kg/m³
    # Kn = λ/r₁ varied by changing gas pressure (thus changing λ and C_s)
    #
    # Table 2's β = K/K_o_base where K_o_base = 8kT/(3η) is the Smoluchowski
    # value with C_s=1. K = K_o_actual * β₂ where K_o_actual includes C_s
    # through diffusion coefficients. So β = (K_o_actual/K_o_base) * β₂ = β₁ * β₂.
    # For equal particles: β₁ = C_s, so β = C_s * β₂.
    sys = DahnekeCoagulationRate()
    compiled = mtkcompile(sys)

    T_val = 298.15  # 25°C
    μ_val = 1.84e-5  # Pa·s
    r_val = 0.22e-6  # 0.22 μm
    ρ_val = 917.0  # kg/m³
    δ_val = 1.0

    # Table 2 data: (Kn, C_s, β_present_theory)
    table2 = [
        (0.1, 1.123, 1.096),
        (0.2, 1.248, 1.214),
        (0.5, 1.653, 1.594),
        (1.0, 2.406, 2.282),
        (2.0, 4.002, 3.661),
        (5.0, 8.907, 7.281),
        (10.0, 17.13, 11.65),
        (20.0, 33.61, 16.54),
        (50.0, 83.04, 21.05),
        (100.0, 165.4, 22.46),
        (200.0, 330.2, 22.95),
        (1000.0, 1649.0, 23.14),
    ]

    k_B = 1.380649e-23
    K_o_smoluchowski = 8 * k_B * T_val / (3 * μ_val)

    for (Kn, C_s, β_expected) in table2
        prob = NonlinearProblem(compiled,
            Dict(compiled.T => T_val, compiled.μ => μ_val,
                compiled.r_1 => r_val, compiled.r_2 => r_val,
                compiled.ρ_1 => ρ_val, compiled.ρ_2 => ρ_val,
                compiled.C_s1 => C_s, compiled.C_s2 => C_s,
                compiled.δ_p => δ_val))
        sol = solve(prob)

        # β = K / K_o_smoluchowski
        # K = K_o_actual * β₂ where K_o_actual = 4π(2r)(2kTC_s/(6πηr)) = 8kTC_s/(3η)
        # So β = C_s * β₂
        β_computed = sol[compiled.K] / K_o_smoluchowski

        @test β_computed ≈ β_expected rtol = 0.02
    end
end

@testitem "DahnekeCoagulationRate δ dependence" setup=[DahnekeSetup] tags=[:dahneke] begin
    # Lower sticking probability should reduce the coagulation rate
    sys = DahnekeCoagulationRate()
    compiled = mtkcompile(sys)

    prob_δ1 = NonlinearProblem(compiled,
        Dict(compiled.r_1 => 1e-7, compiled.r_2 => 1e-7, compiled.δ_p => 1.0))
    sol_δ1 = solve(prob_δ1)

    prob_δ01 = NonlinearProblem(compiled,
        Dict(compiled.r_1 => 1e-7, compiled.r_2 => 1e-7, compiled.δ_p => 0.01))
    sol_δ01 = solve(prob_δ01)

    @test sol_δ01[compiled.K] < sol_δ1[compiled.K]
end

# ============================================================
# DahnekeCapillaryPenetration tests
# ============================================================

@testitem "DahnekeCapillaryPenetration structure" setup=[DahnekeSetup] tags=[:dahneke] begin
    sys = DahnekeCapillaryPenetration()
    @test sys isa System
    @test length(equations(sys)) == 2  # σ_flow, z_dim
end

@testitem "DahnekeCapillaryPenetration σ formula" setup=[DahnekeSetup] tags=[:dahneke] begin
    # Eq. 9.5: σ = 2/(2-γ)
    # γ=0 → σ=1 (plug flow), γ=0.5 → σ=4/3, γ=1.0 → σ=2 (Poiseuille)
    sys = DahnekeCapillaryPenetration()
    compiled = mtkcompile(sys)

    for (γ, σ_expected) in [(0.0, 1.0), (0.5, 4.0 / 3.0), (1.0, 2.0)]
        prob = NonlinearProblem(compiled, Dict(compiled.γ_flow => γ))
        sol = solve(prob)
        @test sol[compiled.σ_flow] ≈ σ_expected rtol = 1e-10
    end
end

@testitem "DahnekeCapillaryPenetration z_dim formula" setup=[DahnekeSetup] tags=[:dahneke] begin
    # z' = DL/(σv̄R²)
    sys = DahnekeCapillaryPenetration()
    compiled = mtkcompile(sys)

    D_val = 1e-10
    L_val = 0.1
    R_val = 1e-4
    v_val = 0.01
    γ_val = 1.0
    σ_val = 2.0  # For γ=1

    prob = NonlinearProblem(compiled,
        Dict(compiled.Dcoeff => D_val, compiled.L_cap => L_val,
            compiled.R_cap => R_val, compiled.v_mean => v_val,
            compiled.γ_flow => γ_val))
    sol = solve(prob)

    z_expected = D_val * L_val / (σ_val * v_val * R_val^2)
    @test sol[compiled.z_dim] ≈ z_expected rtol = 1e-10
end
