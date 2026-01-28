@testitem "ECTracerMethod structure" tags=[:organic_aerosol] begin
    using ModelingToolkit
    using Aerosol

    sys = ECTracerMethod()
    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 4
    @test length(parameters(sys)) == 3  # 2 parameters + 1 constant
end

@testitem "ECTracerMethod equations" tags=[:organic_aerosol] begin
    using ModelingToolkit, NonlinearSolve
    using Aerosol

    sys = ECTracerMethod()
    sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
    ssys = mtkcompile(sys; inputs=[sys_nns.OC, sys_nns.EC])

    # Test with typical values: OC = 10 μg/m³ = 10e-9 kg/m³, EC = 3 μg/m³ = 3e-9 kg/m³
    # Primary OC = 1.7 * 3e-9 + 0.9e-9 = 6.0e-9 kg/m³
    # Secondary OC = 10e-9 - 6.0e-9 = 4.0e-9 kg/m³
    prob = NonlinearProblem(ssys, Dict(
        ssys.OC => 10.0e-9, ssys.EC => 3.0e-9,
        ssys.OC_primary => 1.0e-9, ssys.OC_secondary => 1.0e-9))
    sol = solve(prob)
    @test sol[ssys.OC_primary] ≈ 6.0e-9 rtol=1e-6
    @test sol[ssys.OC_secondary] ≈ 4.0e-9 rtol=1e-6

    # Test with low OC: secondary should be clamped to zero
    prob2 = NonlinearProblem(ssys, Dict(
        ssys.OC => 5.0e-9, ssys.EC => 5.0e-9,
        ssys.OC_primary => 1.0e-9, ssys.OC_secondary => 1.0e-9))
    sol2 = solve(prob2)
    @test sol2[ssys.OC_secondary] ≈ 0.0 atol=1e-18
end

@testitem "NoninteractingSOA structure" tags=[:organic_aerosol] begin
    using ModelingToolkit
    using Aerosol

    sys = NoninteractingSOA()
    @test length(equations(sys)) == 7
    @test length(unknowns(sys)) == 8
    @test length(parameters(sys)) == 8  # 5 parameters + 3 constants
end

@testitem "NoninteractingSOA equations" tags=[:organic_aerosol] begin
    using ModelingToolkit, NonlinearSolve
    using Aerosol

    sys = NoninteractingSOA()
    sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
    ssys = mtkcompile(sys; inputs=[sys_nns.ΔROG])

    # Default parameters in SI
    R = 8.314
    T = 298.0
    p_i = 1.01325e-5
    M_i = 0.180
    M_ROG = 0.150
    a_i = 0.05

    c_eq_expected = p_i * M_i / (R * T)
    threshold_expected = p_i * M_ROG / (a_i * R * T)

    # Test above threshold
    ΔROG_val = 10 * threshold_expected
    c_total_expected = a_i * (M_i / M_ROG) * ΔROG_val
    c_aer_expected = c_total_expected - c_eq_expected

    prob = NonlinearProblem(ssys, Dict(
        ssys.ΔROG => ΔROG_val,
        ssys.c_eq => 1e-10, ssys.c_total => 1e-10,
        ssys.c_aer => 1e-10, ssys.c_gas => 1e-10,
        ssys.ΔROG_threshold => 1e-10,
        ssys.X_p => 0.5, ssys.Y => 0.01))
    sol = solve(prob)
    @test sol[ssys.c_eq] ≈ c_eq_expected rtol=1e-6
    @test sol[ssys.c_total] ≈ c_total_expected rtol=1e-6
    @test sol[ssys.c_aer] ≈ c_aer_expected rtol=1e-6
    @test sol[ssys.c_gas] ≈ c_eq_expected rtol=1e-6
    @test sol[ssys.ΔROG_threshold] ≈ threshold_expected rtol=1e-6
    @test sol[ssys.Y] ≈ c_aer_expected / ΔROG_val rtol=1e-6

    # Test below threshold: all product in gas phase
    ΔROG_below = 0.1 * threshold_expected
    c_total_below = a_i * (M_i / M_ROG) * ΔROG_below
    prob2 = NonlinearProblem(ssys, Dict(
        ssys.ΔROG => ΔROG_below,
        ssys.c_eq => 1e-10, ssys.c_total => 1e-10,
        ssys.c_aer => 1e-10, ssys.c_gas => 1e-10,
        ssys.ΔROG_threshold => 1e-10,
        ssys.X_p => 0.5, ssys.Y => 0.01))
    sol2 = solve(prob2)
    @test sol2[ssys.c_aer] ≈ 0.0 atol=1e-20
    @test sol2[ssys.c_gas] ≈ c_total_below rtol=1e-6
    @test sol2[ssys.X_p] ≈ 0.0 atol=1e-15
end

@testitem "AbsorptivePartitioning structure" tags=[:organic_aerosol] begin
    using ModelingToolkit
    using Aerosol

    sys = AbsorptivePartitioning()
    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 4
    @test length(parameters(sys)) == 10  # 7 parameters + 3 constants
end

@testitem "AbsorptivePartitioning equations" tags=[:organic_aerosol] begin
    using ModelingToolkit, NonlinearSolve
    using Aerosol

    sys = AbsorptivePartitioning()
    sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
    ssys = mtkcompile(sys; inputs=[sys_nns.ΔROG])

    # Default parameters in SI
    R = 8.314
    T = 298.0
    a_i = 0.05
    p_i = 1.01325e-5
    M_i = 0.180
    M_ROG = 0.150
    M_0 = 0.200
    m_0 = 10.0e-9

    X_p_expected = m_0 * R * T / (m_0 * R * T + p_i * M_0)

    ΔROG_val = 50.0e-9
    c_aer_expected = (a_i * R * T / M_ROG) * (M_i * m_0 / (m_0 * R * T + M_0 * p_i)) * ΔROG_val

    prob = NonlinearProblem(ssys, Dict(
        ssys.ΔROG => ΔROG_val,
        ssys.c_aer => 1e-10, ssys.X_p => 0.5, ssys.Y => 0.01))
    sol = solve(prob)
    @test sol[ssys.X_p] ≈ X_p_expected rtol=1e-6
    @test sol[ssys.c_aer] ≈ c_aer_expected rtol=1e-6
    @test sol[ssys.Y] ≈ c_aer_expected / ΔROG_val rtol=1e-6

    # Test linearity: doubling ΔROG should double c_aer
    prob2 = NonlinearProblem(ssys, Dict(
        ssys.ΔROG => 2 * ΔROG_val,
        ssys.c_aer => 1e-10, ssys.X_p => 0.5, ssys.Y => 0.01))
    sol2 = solve(prob2)
    @test sol2[ssys.c_aer] ≈ 2 * c_aer_expected rtol=1e-6
end

@testitem "AbsorptivePartitioning limiting behavior" tags=[:organic_aerosol] begin
    using ModelingToolkit, NonlinearSolve
    using Aerosol

    sys = AbsorptivePartitioning()
    sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
    ssys = mtkcompile(sys; inputs=[sys_nns.ΔROG])

    # When p_i → 0 (very involatile compound), X_p → 1 (everything in aerosol)
    prob_low_p = NonlinearProblem(ssys, Dict(
        ssys.ΔROG => 50.0e-9,
        ssys.c_aer => 1e-10, ssys.X_p => 0.5, ssys.Y => 0.01,
        ssys.p_i => 1e-20))
    sol_low_p = solve(prob_low_p)
    @test sol_low_p[ssys.X_p] ≈ 1.0 atol=1e-6

    # When m_0 → 0, X_p → 0 (no absorbing medium)
    prob_no_oa = NonlinearProblem(ssys, Dict(
        ssys.ΔROG => 50.0e-9,
        ssys.c_aer => 1e-10, ssys.X_p => 0.5, ssys.Y => 0.01,
        ssys.m_0 => 1e-20))
    sol_no_oa = solve(prob_no_oa)
    @test sol_no_oa[ssys.X_p] ≈ 0.0 atol=1e-6
end

@testitem "TwoProductSOA structure" tags=[:organic_aerosol] begin
    using ModelingToolkit
    using Aerosol

    sys = TwoProductSOA()
    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 4
    @test length(parameters(sys)) == 4
end

@testitem "TwoProductSOA equations" tags=[:organic_aerosol] begin
    using ModelingToolkit, NonlinearSolve
    using Aerosol

    sys = TwoProductSOA()
    sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
    ssys = mtkcompile(sys; inputs=[sys_nns.ΔROG])

    # Default: α-pinene/OH parameters from Table 14.12
    a_1 = 0.038
    a_2 = 0.326
    c_sat_1 = 5.8e-9
    c_sat_2 = 250.0e-9

    # Eq. 14.39: threshold = 1 / (a_1/c_sat_1 + a_2/c_sat_2)
    threshold_expected = 1.0 / (a_1 / c_sat_1 + a_2 / c_sat_2)

    prob = NonlinearProblem(ssys, Dict(
        ssys.ΔROG => 100.0e-9,
        ssys.c_aer => 10.0e-9,
        ssys.Y => 0.1, ssys.ΔROG_threshold => 1e-9))
    sol = solve(prob)
    @test sol[ssys.ΔROG_threshold] ≈ threshold_expected rtol=1e-6

    # Verify Eq. 14.43 at the solved c_aer value
    c_aer_val = sol[ssys.c_aer]
    Y_expected = c_aer_val * (a_1 / (c_aer_val + c_sat_1) + a_2 / (c_aer_val + c_sat_2))
    @test sol[ssys.Y] ≈ Y_expected rtol=1e-6
end

@testitem "LangmuirAdsorption structure" tags=[:organic_aerosol] begin
    using ModelingToolkit
    using Aerosol

    sys = LangmuirAdsorption()
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 2
    @test length(parameters(sys)) == 2
end

@testitem "LangmuirAdsorption equations" tags=[:organic_aerosol] begin
    using ModelingToolkit, NonlinearSolve
    using Aerosol

    sys = LangmuirAdsorption()
    sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
    ssys = mtkcompile(sys; inputs=[sys_nns.p])

    V_m = 1.0
    b = 1.0

    # Test at p = 1 Pa: V = V_m * b * p / (1 + b*p) = 1*1*1/(1+1) = 0.5
    prob = NonlinearProblem(ssys, Dict(ssys.p => 1.0, ssys.V => 0.5))
    sol = solve(prob)
    @test sol[ssys.V] ≈ 0.5 rtol=1e-6

    # At very high pressure: V → V_m (saturation)
    prob_high = NonlinearProblem(ssys, Dict(ssys.p => 1e6, ssys.V => 0.5))
    sol_high = solve(prob_high)
    @test sol_high[ssys.V] ≈ V_m rtol=1e-3

    # At very low pressure: V ≈ V_m * b * p (Henry's law limit)
    p_low = 1e-6
    prob_low = NonlinearProblem(ssys, Dict(ssys.p => p_low, ssys.V => 0.5))
    sol_low = solve(prob_low)
    @test sol_low[ssys.V] ≈ V_m * b * p_low rtol=1e-3
end

@testitem "BETAdsorption structure" tags=[:organic_aerosol] begin
    using ModelingToolkit
    using Aerosol

    sys = BETAdsorption()
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 2
    @test length(parameters(sys)) == 2
end

@testitem "BETAdsorption equations" tags=[:organic_aerosol] begin
    using ModelingToolkit, NonlinearSolve
    using Aerosol

    sys = BETAdsorption()
    sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
    ssys = mtkcompile(sys; inputs=[sys_nns.S])

    V_m = 1.0
    c_BET = 10.0

    S_val = 0.5
    V_expected = V_m * c_BET * S_val / ((1 - S_val) * (1 + (c_BET - 1) * S_val))

    prob = NonlinearProblem(ssys, Dict(ssys.S => S_val, ssys.V => 1.0))
    sol = solve(prob)
    @test sol[ssys.V] ≈ V_expected rtol=1e-6

    # As S → 0 (low humidity), V → V_m * c_BET * S (linear in S)
    S_low = 0.001
    V_low_expected = V_m * c_BET * S_low / ((1 - S_low) * (1 + (c_BET - 1) * S_low))
    prob_low = NonlinearProblem(ssys, Dict(ssys.S => S_low, ssys.V => 0.01))
    sol_low = solve(prob_low)
    @test sol_low[ssys.V] ≈ V_low_expected rtol=1e-6
    # At low S, approximately V_m * c_BET * S
    @test sol_low[ssys.V] ≈ V_m * c_BET * S_low rtol=1e-2
end

@testitem "FHHAdsorption structure" tags=[:organic_aerosol] begin
    using ModelingToolkit
    using Aerosol

    sys = FHHAdsorption()
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 2
    @test length(parameters(sys)) == 3
end

@testitem "FHHAdsorption equations" tags=[:organic_aerosol] begin
    using ModelingToolkit, NonlinearSolve
    using Aerosol

    sys = FHHAdsorption()
    sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
    ssys = mtkcompile(sys; inputs=[sys_nns.S])

    V_m = 1.0
    A_FHH = 1.0
    B_FHH = 1.0

    # Eq. 14.46: V = V_m * (A / ln(1/S))^(1/B)
    S_val = 0.5
    V_expected = V_m * (A_FHH / log(1 / S_val))^(1 / B_FHH)

    prob = NonlinearProblem(ssys, Dict(ssys.S => S_val, ssys.V => 1.0))
    sol = solve(prob)
    @test sol[ssys.V] ≈ V_expected rtol=1e-6
    @test sol[ssys.V] ≈ 1.0 / log(2.0) rtol=1e-6

    # V increases as S approaches 1
    S_high = 0.99
    prob_high = NonlinearProblem(ssys, Dict(ssys.S => S_high, ssys.V => 10.0))
    sol_high = solve(prob_high)
    @test sol_high[ssys.V] > sol[ssys.V]

    # Test with different A,B parameters
    prob_b2 = NonlinearProblem(ssys, Dict(
        ssys.S => 0.5, ssys.V => 1.0,
        ssys.A_FHH => 2.0, ssys.B_FHH => 2.0))
    sol_b2 = solve(prob_b2)
    V_expected_b2 = V_m * (2.0 / log(2.0))^0.5
    @test sol_b2[ssys.V] ≈ V_expected_b2 rtol=1e-6
end

@testitem "BET reduces to Langmuir" tags=[:organic_aerosol] begin
    using ModelingToolkit, NonlinearSolve
    using Aerosol

    bet = BETAdsorption()
    bet_nns = ModelingToolkit.toggle_namespacing(bet, false)
    sbet = mtkcompile(bet; inputs=[bet_nns.S])

    lang = LangmuirAdsorption()
    lang_nns = ModelingToolkit.toggle_namespacing(lang, false)
    slang = mtkcompile(lang; inputs=[lang_nns.p])

    S_low = 0.001
    c_BET = 10.0

    bet_prob = NonlinearProblem(sbet, Dict(sbet.S => S_low, sbet.V => 0.01))
    bet_sol = solve(bet_prob)

    # Equivalent Langmuir: b*p = c_BET * S
    lang_prob = NonlinearProblem(slang, Dict(slang.p => c_BET * S_low, slang.V => 0.01))
    lang_sol = solve(lang_prob)

    @test bet_sol[sbet.V] ≈ lang_sol[slang.V] rtol=1e-2
end
