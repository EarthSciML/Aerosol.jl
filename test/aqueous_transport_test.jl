@testsnippet AqueousTransportSetup begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEq
    using NonlinearSolve
    using Aerosol
end

@testitem "AqueousDiffusionReaction structure" setup=[AqueousTransportSetup] tags=[:aqueous_transport] begin
    sys = AqueousDiffusionReaction()
    @test sys isa System
    @test length(equations(sys)) == 2
end

@testitem "AqueousDiffusionReaction Q limiting values" setup=[AqueousTransportSetup] tags=[:aqueous_transport] begin
    sys = AqueousDiffusionReaction()
    compiled = mtkcompile(sys)

    # For small q (slow reaction), Q → 1
    R_p = 1.0e-5
    D_aq = 1.0e-9
    k_slow = 1.0  # s⁻¹ (slow reaction)
    q_slow = R_p * sqrt(k_slow / D_aq)

    # For large q (fast reaction), Q → 3/q (diffusion-limited)
    k_fast = 1.0e6  # s⁻¹ (very fast reaction)
    q_fast = R_p * sqrt(k_fast / D_aq)

    prob_slow = NonlinearProblem(compiled, Dict(compiled.R_p => R_p, compiled.D_aq => D_aq, compiled.k_rxn => k_slow))
    sol_slow = solve(prob_slow)

    prob_fast = NonlinearProblem(compiled, Dict(compiled.R_p => R_p, compiled.D_aq => D_aq, compiled.k_rxn => k_fast))
    sol_fast = solve(prob_fast)

    q_calc_slow = sol_slow[compiled.q]
    Q_slow = sol_slow[compiled.Q]

    q_calc_fast = sol_fast[compiled.q]
    Q_fast = sol_fast[compiled.Q]

    # Verify q calculation
    @test q_calc_slow ≈ q_slow rtol=1e-6
    @test q_calc_fast ≈ q_fast rtol=1e-6

    # For slow reaction (small q), Q should be close to 1
    # Note: The formula Q = 3(coth(q)/q - 1/q²) has a numerical limit issue for very small q
    # For q ≈ 0.316, Q should be close to 1
    @test Q_slow > 0.8

    # For fast reaction (large q), Q ≈ 3/q
    expected_Q_fast = 3 / q_fast
    @test Q_fast ≈ expected_Q_fast rtol=0.1
end

@testitem "MassTransportLimitation structure" setup=[AqueousTransportSetup] tags=[:aqueous_transport] begin
    sys = MassTransportLimitation()
    @test sys isa System
    @test length(equations(sys)) == 3
end

@testitem "MassTransportLimitation scaling" setup=[AqueousTransportSetup] tags=[:aqueous_transport] begin
    sys = MassTransportLimitation()
    compiled = mtkcompile(sys)

    R_p = 1.0e-5  # 10 μm
    D_g = 2.0e-5
    D_aq = 1.0e-9
    T = 298.15
    α = 1.0
    M_A = 0.029

    prob = NonlinearProblem(compiled, Dict(compiled.R_p => R_p, compiled.D_g => D_g, compiled.D_aq => D_aq,
         compiled.T => T, compiled.α => α, compiled.M_A => M_A))
    sol = solve(prob)

    k1H_gas = sol[compiled.k1H_gas_limit]
    k1_aq = sol[compiled.k1_aq_limit]
    k1H_int = sol[compiled.k1H_interface_limit]

    # All limits should be positive
    @test k1H_gas > 0
    @test k1_aq > 0
    @test k1H_int > 0

    # Scaling with R_p: gas limit ∝ R_p⁻², aqueous limit ∝ R_p⁻², interface ∝ R_p⁻¹
    # Test by comparing two droplet sizes
    R_p2 = 2.0e-5

    prob2 = NonlinearProblem(compiled, Dict(compiled.R_p => R_p2, compiled.D_g => D_g, compiled.D_aq => D_aq,
         compiled.T => T, compiled.α => α, compiled.M_A => M_A))
    sol2 = solve(prob2)

    k1H_gas_2 = sol2[compiled.k1H_gas_limit]
    k1_aq_2 = sol2[compiled.k1_aq_limit]
    k1H_int_2 = sol2[compiled.k1H_interface_limit]

    # Gas and aqueous limits should scale as (R_p/R_p2)²
    @test k1H_gas / k1H_gas_2 ≈ (R_p2 / R_p)^2 rtol=1e-6
    @test k1_aq / k1_aq_2 ≈ (R_p2 / R_p)^2 rtol=1e-6

    # Interface limit should scale as R_p/R_p2
    @test k1H_int / k1H_int_2 ≈ R_p2 / R_p rtol=1e-6
end

@testitem "DropletMassBalance structure" setup=[AqueousTransportSetup] tags=[:aqueous_transport] begin
    sys = DropletMassBalance()
    @test sys isa System
    @test length(equations(sys)) == 2
end

@testitem "DropletMassBalance dynamics" setup=[AqueousTransportSetup] tags=[:aqueous_transport] begin
    # Test that the mass balance equations produce sensible dynamics:
    # - When gas phase has higher concentration than equilibrium with aqueous,
    #   gas should decrease and aqueous should increase (uptake)
    sys = DropletMassBalance()
    compiled = mtkcompile(sys)

    k_mt = 10.0
    w_L = 1.0e-3
    H_star = 1.0e2  # mol/(m³·Pa)
    T = 298.15
    R = 8.314

    # Initial conditions: gas phase has pressure, aqueous phase empty
    p0 = 100.0  # Pa
    C_aq0 = 0.0

    prob = ODEProblem(compiled,
        merge(Dict(compiled.p => p0, compiled.C_aq => C_aq0),
              Dict(compiled.k_mt => k_mt, compiled.w_L => w_L,
                   compiled.H_star => H_star, compiled.T => T,
                   compiled.Q => 1.0, compiled.R_aq => 0.0)),
        (0.0, 1.0))
    sol = solve(prob)

    p_final = sol[compiled.p][end]
    C_aq_final = sol[compiled.C_aq][end]

    # Gas phase pressure should decrease (uptake into aqueous)
    @test p_final < p0

    # Aqueous concentration should increase (uptake from gas)
    @test C_aq_final > C_aq0

    # The system should be approaching equilibrium: p ≈ C_aq / H*
    # After some time, the ratio C_aq / (H* × p) should be closer to 1 than initially
    # Initial ratio is 0/anything = 0
    # At equilibrium, ratio = 1
    # So ratio should be increasing toward 1
    ratio = C_aq_final / (H_star * p_final)
    @test ratio > 0  # Has moved toward equilibrium
    @test ratio < 1  # Hasn't overshot
end

@testitem "DropletMassBalance mass conservation" setup=[AqueousTransportSetup] tags=[:aqueous_transport] begin
    # Test mass conservation when there is no reaction
    sys = DropletMassBalance()
    compiled = mtkcompile(sys)

    k_mt = 1.0
    w_L = 1.0e-6
    H_star = 1.0e5
    T = 298.15
    R = 8.314

    # Initial conditions
    p0 = 100.0  # Pa
    C_aq0 = 0.0

    prob = ODEProblem(compiled,
        merge(Dict(compiled.p => p0, compiled.C_aq => C_aq0),
              Dict(compiled.k_mt => k_mt, compiled.w_L => w_L,
                   compiled.H_star => H_star, compiled.T => T,
                   compiled.Q => 1.0, compiled.R_aq => 0.0)),
        (0.0, 100.0))
    sol = solve(prob)

    # Total moles should be conserved
    # n_gas = p/(RT), n_aq = C_aq × w_L (per unit total volume)
    n_initial = p0 / (R * T) + C_aq0 * w_L

    p_final = sol[compiled.p][end]
    C_aq_final = sol[compiled.C_aq][end]
    n_final = p_final / (R * T) + C_aq_final * w_L

    @test n_final ≈ n_initial rtol=0.01
end
