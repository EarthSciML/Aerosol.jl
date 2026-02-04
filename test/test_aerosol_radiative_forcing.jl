@testsnippet AerosolRadiativeForcingSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D, mtkcompile, unknowns, equations, parameters
    using DynamicQuantities
    using SteadyStateDiffEq
    using Aerosol

    """
    Solve a ModelingToolkit algebraic system for steady-state values.
    Uses NonlinearSolve via SteadyStateDiffEq for systems without dynamics.
    """
    function solve_system(sys; param_vals=Dict())
        csys = mtkcompile(sys)
        # Build parameter map using the compiled system's parameter symbols
        p = Dict{Any, Any}()
        for (name, val) in param_vals
            param_sym = getproperty(csys, name)
            p[param_sym] = val
        end
        # Use the new API: SteadyStateProblem(sys, merge(Dict(), Dict(p)))
        prob = SteadyStateProblem(csys, p)
        sol = solve(prob)
        return sol, csys
    end
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
    # Figure 24.2 (Seinfeld & Pandis 2006, p. 1061) shows curves for different β values

    # Test case: R_s = 0.15, β = 0.29 (Table 24.2 values, p. 1069)
    # Expected ω_crit ≈ 0.6 from Section 24.2 (p. 1061)
    sys = CriticalSingleScatteringAlbedo()
    sol, csys = solve_system(sys, param_vals=Dict(:R_s => 0.15, :β => 0.29))

    # Verify the system computes expected value
    ω_crit_computed = sol[csys.ω_crit][1]
    @test ω_crit_computed ≈ 0.586 rtol=0.01  # Approximately 0.6

    # Test limiting case from Figure 24.2: At R_s → 0: ω_crit → 0
    sol_low, csys = solve_system(CriticalSingleScatteringAlbedo(), param_vals=Dict(:R_s => 0.01, :β => 0.29))
    @test sol_low[csys.ω_crit][1] < 0.1

    # Test limiting case: At R_s → 1: ω_crit → 1
    sol_high, csys = solve_system(CriticalSingleScatteringAlbedo(), param_vals=Dict(:R_s => 0.99, :β => 0.29))
    @test sol_high[csys.ω_crit][1] > 0.99

    # Verify ω_crit increases with increasing R_s
    R_s_vals = [0.1, 0.3, 0.5, 0.7, 0.9]
    ω_crit_vals = Float64[]
    for R_s in R_s_vals
        sol_i, csys_i = solve_system(CriticalSingleScatteringAlbedo(), param_vals=Dict(:R_s => R_s, :β => 0.29))
        push!(ω_crit_vals, sol_i[csys_i.ω_crit][1])
    end
    @test issorted(ω_crit_vals)

    # Verify ω_crit decreases with increasing β (at fixed R_s)
    β_vals = [0.1, 0.2, 0.3, 0.4]
    ω_crit_vs_beta = Float64[]
    for β in β_vals
        sol_i, csys_i = solve_system(CriticalSingleScatteringAlbedo(), param_vals=Dict(:R_s => 0.5, :β => β))
        push!(ω_crit_vs_beta, sol_i[csys_i.ω_crit][1])
    end
    @test issorted(ω_crit_vs_beta, rev=true)
end

@testitem "Cloud Albedo - Figure 24.16 Validation" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test Eq. 24.37-24.38 (Seinfeld & Pandis 2006, p. 1081)
    # R_c = τ_c / (τ_c + γ) where γ = 2/(√3*(1-g)) ≈ 7.7 for g = 0.85

    # First verify γ is computed correctly by the system
    sol, csys = solve_system(CloudAlbedo(), param_vals=Dict(:τ_c => 10.0, :g => 0.85))
    γ_computed = sol[csys.γ][1]
    @test γ_computed ≈ 7.698 rtol=0.01

    # Test key points from Figure 24.16 (p. 1082)
    # At τ_c ≈ 0: R_c → 0
    sol_zero, csys = solve_system(CloudAlbedo(), param_vals=Dict(:τ_c => 0.01, :g => 0.85))
    @test sol_zero[csys.R_c][1] < 0.01

    # At τ_c = γ: R_c = 0.5 (half-albedo point)
    sol_half, csys = solve_system(CloudAlbedo(), param_vals=Dict(:τ_c => γ_computed, :g => 0.85))
    @test sol_half[csys.R_c][1] ≈ 0.5 rtol=0.01

    # At τ_c >> γ: R_c → 1
    sol_large, csys = solve_system(CloudAlbedo(), param_vals=Dict(:τ_c => 100.0, :g => 0.85))
    @test sol_large[csys.R_c][1] > 0.9

    # Verify R_c increases monotonically with τ_c
    τ_c_vals = [1.0, 5.0, 10.0, 20.0, 50.0]
    R_c_vals = Float64[]
    for τ_c in τ_c_vals
        sol_i, csys_i = solve_system(CloudAlbedo(), param_vals=Dict(:τ_c => τ_c, :g => 0.85))
        push!(R_c_vals, sol_i[csys_i.R_c][1])
    end
    @test issorted(R_c_vals)
end

@testitem "Twomey Susceptibility - Eq. 24.40-24.41 Validation" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test Eq. 24.40-24.41 (Seinfeld & Pandis 2006, p. 1083)
    # dR_c/d(ln N) = R_c*(1-R_c)/3 (Twomey susceptibility S)

    # Maximum susceptibility occurs at R_c = 0.5
    sol, csys = solve_system(CloudAlbedoSensitivity(), param_vals=Dict(:R_c => 0.5, :N => 100.0e6))
    S_max = sol[csys.S][1]
    @test S_max ≈ 0.0833 rtol=0.01  # 1/12 ≈ 0.0833

    # Also verify dRc_dN is computed correctly
    dRc_dN_computed = sol[csys.dRc_dN][1]
    @test dRc_dN_computed ≈ (0.5 * 0.5 / 3) / 100.0e6 rtol=1e-6

    # Susceptibility at R_c = 0.5 is maximum - verify by testing multiple R_c values
    R_c_vals = [0.1, 0.3, 0.5, 0.7, 0.9]
    S_vals = Float64[]
    for R_c in R_c_vals
        sol_i, csys_i = solve_system(CloudAlbedoSensitivity(), param_vals=Dict(:R_c => R_c, :N => 100.0e6))
        push!(S_vals, sol_i[csys_i.S][1])
    end

    # Find index of maximum
    max_idx = argmax(S_vals)
    @test R_c_vals[max_idx] ≈ 0.5

    # Susceptibility is symmetric around R_c = 0.5
    @test S_vals[1] ≈ S_vals[5] rtol=0.01  # S(0.1) ≈ S(0.9)
    @test S_vals[2] ≈ S_vals[4] rtol=0.01  # S(0.3) ≈ S(0.7)
end

@testitem "Indirect Forcing - Section 24.8.2 Validation" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test Eq. 24.42-24.43 (Seinfeld & Pandis 2006, p. 1083)
    # ΔF_c = -F_0 * A_c * T_a² * ΔR_c

    # Parameters from Table 24.2 (p. 1069)
    # 30% increase in N: Δln(N) = ln(1.3) ≈ 0.262
    Δln_N = log(1.3)
    @test Δln_N ≈ 0.262 rtol=0.01

    # Test with the implementation
    sol, csys = solve_system(IndirectAerosolForcing(), param_vals=Dict(
        :F_0 => 1370.0,
        :A_c => 0.6,
        :T_a => 0.76,
        :R_c => 0.5,
        :Δln_N => Δln_N
    ))

    # Verify ΔR_c computed by system
    ΔR_c_computed = sol[csys.ΔR_c][1]
    @test ΔR_c_computed ≈ 0.0219 rtol=0.05

    # Verify ΔF_c computed by system
    ΔF_c_computed = sol[csys.ΔF_c][1]

    # This gives the cloud forcing sensitivity for a 30% change in N globally
    # The magnitude is consistent with ERBE observations that cloud shortwave
    # forcing is ~50 W/m² (Section 24.8.2, p. 1084)
    @test ΔF_c_computed < 0  # Negative = cooling (increased albedo reflects more sunlight)
    @test abs(ΔF_c_computed) > 5.0  # Significant forcing magnitude
    @test abs(ΔF_c_computed) < 20.0  # But bounded to realistic range

    # Verify formula consistency using computed values
    F_0, A_c, T_a = 1370.0, 0.6, 0.76
    ΔF_c_expected = -F_0 * A_c * T_a^2 * ΔR_c_computed
    @test ΔF_c_computed ≈ ΔF_c_expected rtol=1e-6
end

@testitem "Aerosol Layer Forcing - Limiting Cases" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test limiting behavior from Section 24.1 (Seinfeld & Pandis 2006, p. 1057-1059)
    # Using the actual AerosolLayerRadiativeForcing system

    # Case 1: τ = 0 (no aerosol) - reflectance should match surface reflectance
    # When τ → 0, r_aer → 0 and t_aer → 1, so R_as → R_s
    sol_zero, csys = solve_system(AerosolLayerRadiativeForcing(), param_vals=Dict(
        :τ => 1e-10, :ω => 0.95, :β => 0.29, :R_s => 0.15, :A_c => 0.6, :T_a => 0.76, :F_0 => 1370.0
    ))

    r_zero = sol_zero[csys.r_aer][1]
    t_zero = sol_zero[csys.t_aer][1]
    ΔR_p_zero = sol_zero[csys.ΔR_p][1]

    @test r_zero ≈ 0.0 atol=1e-8
    @test t_zero ≈ 1.0 atol=1e-8
    @test ΔR_p_zero ≈ 0.0 atol=1e-8  # No aerosol = no change in albedo

    # Case 2: ω = 1 (nonabsorbing aerosol)
    # For nonabsorbing aerosol: r + t should equal 1 (energy conservation)
    sol_nonabs, csys = solve_system(AerosolLayerRadiativeForcing(), param_vals=Dict(
        :τ => 0.1, :ω => 1.0, :β => 0.29, :R_s => 0.15, :A_c => 0.6, :T_a => 0.76, :F_0 => 1370.0
    ))

    r_nonabs = sol_nonabs[csys.r_aer][1]
    t_nonabs = sol_nonabs[csys.t_aer][1]
    @test r_nonabs + t_nonabs ≈ 1.0 rtol=0.01

    # Case 3: τ << 1 approximation (Eq. 24.11, p. 1060)
    # r ≈ τωβ for small τ
    τ_small = 0.01
    sol_small, csys = solve_system(AerosolLayerRadiativeForcing(), param_vals=Dict(
        :τ => τ_small, :ω => 0.95, :β => 0.29, :R_s => 0.15, :A_c => 0.6, :T_a => 0.76, :F_0 => 1370.0
    ))

    r_small = sol_small[csys.r_aer][1]
    r_approx = τ_small * 0.95 * 0.29  # Linear approximation for small τ
    @test r_small ≈ r_approx rtol=0.05  # Good approximation for small τ
end

@testitem "Cloud Optical Depth - Eq. 24.36" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test Eq. 24.36 (Seinfeld & Pandis 2006, p. 1081):
    # τ_c = h * (9π*L²*N / 2ρ_w²)^(1/3)

    # Typical values from Figure 24.17 caption (p. 1082)
    sol, csys = solve_system(CloudOpticalDepth(), param_vals=Dict(
        :h => 500.0, :L => 0.3e-3, :N => 100.0e6
    ))
    τ_c = sol[csys.τ_c][1]

    # Cloud optical depth should be positive and reasonable
    @test τ_c > 0
    @test τ_c < 100  # Reasonable upper bound

    # τ_c should increase with N (more droplets = more optical depth)
    sol_high_N, csys = solve_system(CloudOpticalDepth(), param_vals=Dict(
        :h => 500.0, :L => 0.3e-3, :N => 200.0e6
    ))
    @test sol_high_N[csys.τ_c][1] > τ_c

    # τ_c should increase with L
    sol_high_L, csys = solve_system(CloudOpticalDepth(), param_vals=Dict(
        :h => 500.0, :L => 0.5e-3, :N => 100.0e6
    ))
    @test sol_high_L[csys.τ_c][1] > τ_c

    # τ_c should scale linearly with h
    sol_double_h, csys = solve_system(CloudOpticalDepth(), param_vals=Dict(
        :h => 1000.0, :L => 0.3e-3, :N => 100.0e6
    ))
    @test sol_double_h[csys.τ_c][1] ≈ 2 * τ_c rtol=0.01
end

@testitem "Direct Forcing - Cooling vs Heating" setup=[AerosolRadiativeForcingSetup] tags=[:aerosol_radiative] begin
    # Test cooling vs heating boundary (Seinfeld & Pandis 2006, Section 24.2, p. 1060-1061)
    # Aerosols with ω > ω_crit cause cooling, ω < ω_crit cause heating

    # First get ω_crit for typical surface
    sol_crit, csys_crit = solve_system(CriticalSingleScatteringAlbedo(), param_vals=Dict(:R_s => 0.15, :β => 0.29))
    ω_crit = sol_crit[csys_crit.ω_crit][1]

    # Scattering aerosol (ω > ω_crit) - should cool
    ω_scatter = 0.95
    @test ω_scatter > ω_crit

    sol_scatter, csys = solve_system(AerosolLayerRadiativeForcing(), param_vals=Dict(
        :τ => 0.1, :ω => ω_scatter, :β => 0.29, :R_s => 0.15, :A_c => 0.6, :T_a => 0.76, :F_0 => 1370.0
    ))
    ΔR_p_scatter = sol_scatter[csys.ΔR_p][1]
    @test ΔR_p_scatter > 0  # Increased albedo = cooling

    # Absorbing aerosol (ω < ω_crit) - should heat
    ω_absorb = 0.4
    @test ω_absorb < ω_crit

    sol_absorb, csys = solve_system(AerosolLayerRadiativeForcing(), param_vals=Dict(
        :τ => 0.1, :ω => ω_absorb, :β => 0.29, :R_s => 0.15, :A_c => 0.6, :T_a => 0.76, :F_0 => 1370.0
    ))
    ΔR_p_absorb = sol_absorb[csys.ΔR_p][1]
    @test ΔR_p_absorb < 0  # Decreased albedo = heating
end
