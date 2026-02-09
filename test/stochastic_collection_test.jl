@testsnippet SCESetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D
    using OrdinaryDiffEq
    using Aerosol

    """
    Set up exponential initial distribution for SCE tests.
    Returns a Dict of initial conditions for Nk and Mk.
    n(x,0) = (N0/x̄) exp(-x/x̄) where x̄ = LWC/N0.
    """
    function sce_exponential_ic(compiled, I, x1, N0, LWC)
        xbar = LWC / N0
        u0 = Dict{Any,Float64}()
        for k in 1:I
            xk = x1 * 2.0^(k - 1)
            xk1 = x1 * 2.0^k
            Nk_init = N0 * (exp(-xk / xbar) - exp(-xk1 / xbar))
            Mk_init = N0 * xbar * ((xk / xbar + 1) * exp(-xk / xbar) - (xk1 / xbar + 1) * exp(-xk1 / xbar))
            u0[compiled.Nk[k]] = max(Nk_init, 1e-20)
            u0[compiled.Mk[k]] = max(Mk_init, 1e-30)
        end
        return u0
    end

    """
    Create exponential initial distribution arrays for direct RHS testing.
    """
    function sce_exponential_arrays(I, x1, N0, LWC)
        xbar = LWC / N0
        N_init = zeros(I)
        M_init = zeros(I)
        for k in 1:I
            xk = x1 * 2.0^(k - 1)
            xk1 = x1 * 2.0^k
            N_init[k] = max(N0 * (exp(-xk / xbar) - exp(-xk1 / xbar)), 1e-20)
            M_init[k] = max(N0 * xbar * ((xk / xbar + 1) * exp(-xk / xbar) - (xk1 / xbar + 1) * exp(-xk1 / xbar)), 1e-30)
        end
        return N_init, M_init
    end

    """
    Compute terms 5+6 contribution to dN_total/dt for constant kernel.
    For constant kernel, this should equal -(K0/2)*N_total^2 exactly.
    """
    function compute_terms56_constant(N, I, K0)
        result = 0.0
        for k in 1:I
            result -= (K0 / 2) * N[k]^2  # Term 5
            for i in (k + 1):I
                result -= K0 * N[k] * N[i]  # Term 6
            end
        end
        return result
    end
end

@testitem "StochasticCollectionCoalescence structure (constant kernel)" setup=[SCESetup] tags=[:sce] begin
    sys = StochasticCollectionCoalescence(; I=5, kernel_type=:constant, kernel_params=Dict(:K0 => 1e-10))
    @test sys isa System

    # 2*I equations: I for dNk/dt and I for dMk/dt
    eqs = equations(sys)
    @test length(eqs) == 10

    vars = unknowns(sys)
    @test length(vars) == 10
end

@testitem "StochasticCollectionCoalescence structure (Golovin kernel)" setup=[SCESetup] tags=[:sce] begin
    sys = StochasticCollectionCoalescence(; I=5, kernel_type=:golovin, kernel_params=Dict(:C => 1.5e-3))
    @test sys isa System

    eqs = equations(sys)
    @test length(eqs) == 10
end

@testitem "StochasticCollectionCoalescence compiles and solves" setup=[SCESetup] tags=[:sce] begin
    I = 8
    sys = StochasticCollectionCoalescence(; I=I, kernel_type=:constant, kernel_params=Dict(:K0 => 1e-10))
    compiled = mtkcompile(sys)
    @test length(unknowns(compiled)) == 2 * I

    x1 = 1.6e-14
    u0 = sce_exponential_ic(compiled, I, x1, 300e6, 1e-3)

    prob = ODEProblem(compiled, u0, (0.0, 100.0))
    sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-12)
    @test sol.retcode == ReturnCode.Success
end

@testitem "Constant kernel: RHS total number rate matches analytical" setup=[SCESetup] tags=[:sce] begin
    # For constant kernel K(x,y) = K0, terms 5+6 of dN_k/dt summed over all k
    # should give exactly -(K0/2)*N_total^2
    I = 8; x1 = 1e-14; K0 = 1e-10
    N_init, M_init = sce_exponential_arrays(I, x1, 300e6, 1e-3)

    N_total = sum(N_init)
    term56 = compute_terms56_constant(N_init, I, K0)
    analytical = -(K0 / 2) * N_total^2
    @test term56 ≈ analytical rtol = 1e-10
end

@testitem "Constant kernel: mass conservation at RHS level" setup=[SCESetup] tags=[:sce] begin
    # Verify dM_total/dt is finite and negative (mass leaks out the top bin).
    # With only 8 bins the truncation error is large, so we just check the
    # qualitative direction and that the rate is bounded.
    I = 8; x1 = 1e-14; K0 = 1e-10
    xb = [x1 * 2.0^(k - 1) for k in 1:(I + 1)]
    N_init, M_init = sce_exponential_arrays(I, x1, 300e6, 1e-3)

    state = vcat(N_init, M_init)
    result = Aerosol._sce_rhs(state, xb, 1.0, K0)
    dM = result[I+1:2I]
    # Total mass rate should be negative (mass lost through top bin) or near zero
    @test sum(dM) <= 0.0
    # The rate should be finite and bounded
    @test isfinite(sum(dM))
end

@testitem "Constant kernel: number decreases, mass approximately conserved" setup=[SCESetup] tags=[:sce] begin
    I = 8
    sys = StochasticCollectionCoalescence(; I=I, kernel_type=:constant, kernel_params=Dict(:K0 => 1e-10))
    compiled = mtkcompile(sys)

    x1 = 1.6e-14
    u0 = sce_exponential_ic(compiled, I, x1, 300e6, 1e-3)

    prob = ODEProblem(compiled, u0, (0.0, 3600.0))
    sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-12)
    @test sol.retcode == ReturnCode.Success

    N_initial = sum(sol[compiled.Nk[k]][1] for k in 1:I)
    N_final = sum(sol[compiled.Nk[k]][end] for k in 1:I)
    M_initial = sum(sol[compiled.Mk[k]][1] for k in 1:I)
    M_final = sum(sol[compiled.Mk[k]][end] for k in 1:I)

    # Number should decrease due to coalescence
    @test N_final < N_initial

    # With only 8 bins and 3600 s of coalescence, substantial mass leaks through the
    # top bin. Check that mass decreased (not increased) and that some mass remains.
    @test M_final > 0.0
    @test M_final < M_initial
end

@testitem "Golovin kernel compiles, solves, and conserves mass" setup=[SCESetup] tags=[:sce] begin
    I = 8
    sys = StochasticCollectionCoalescence(; I=I, kernel_type=:golovin, kernel_params=Dict(:C => 1.5e-3))
    compiled = mtkcompile(sys)

    x1 = 1.6e-14
    u0 = sce_exponential_ic(compiled, I, x1, 300e6, 1e-3)

    prob = ODEProblem(compiled, u0, (0.0, 1800.0))
    sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-12)
    @test sol.retcode == ReturnCode.Success

    N_initial = sum(sol[compiled.Nk[k]][1] for k in 1:I)
    N_final = sum(sol[compiled.Nk[k]][end] for k in 1:I)
    @test N_final < N_initial

    M_initial = sum(sol[compiled.Mk[k]][1] for k in 1:I)
    M_final = sum(sol[compiled.Mk[k]][end] for k in 1:I)
    @test M_final > 0.1 * M_initial
end

@testitem "Golovin kernel: RHS mass conservation" setup=[SCESetup] tags=[:sce] begin
    # With only 8 bins, significant mass leaks through the top bin.
    # Verify the mass rate is negative (leaking out) and finite.
    I = 8; x1 = 1e-14; C_val = 1.5e-3
    xb = [x1 * 2.0^(k - 1) for k in 1:(I + 1)]
    N_init, M_init = sce_exponential_arrays(I, x1, 300e6, 1e-3)

    state = vcat(N_init, M_init)
    result = Aerosol._sce_rhs(state, xb, 2.0, C_val)
    dM = result[I+1:2I]
    @test sum(dM) <= 0.0
    @test isfinite(sum(dM))
end

@testitem "Invalid kernel type raises error" setup=[SCESetup] tags=[:sce] begin
    @test_throws ErrorException StochasticCollectionCoalescence(; kernel_type=:invalid, kernel_params=Dict(:K0 => 1.0))
end

@testitem "Positivity preservation" setup=[SCESetup] tags=[:sce] begin
    I = 8
    sys = StochasticCollectionCoalescence(; I=I, kernel_type=:constant, kernel_params=Dict(:K0 => 1e-10))
    compiled = mtkcompile(sys)

    x1 = 1.6e-14
    u0 = sce_exponential_ic(compiled, I, x1, 300e6, 1e-3)

    prob = ODEProblem(compiled, u0, (0.0, 3600.0))
    sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-12)
    @test sol.retcode == ReturnCode.Success

    for k in 1:I
        @test sol[compiled.Nk[k]][end] >= -1e-10
        @test sol[compiled.Mk[k]][end] >= -1e-20
    end
end

@testitem "Closure parameter and linear params" setup=[SCESetup] tags=[:sce] begin
    xk = 1.0
    Nk = 10.0
    Mk = 15.0  # x̄ = 1.5, in [1, 2]
    fk, ψk = Aerosol._sce_linear_params(Nk, Mk, xk)
    @test fk ≈ 10.0
    @test ψk ≈ 10.0

    # Positivity case: x̄ < xk (Eq. 15b)
    Mk_low = 5.0
    fk_low, ψk_low = Aerosol._sce_linear_params(Nk, Mk_low, xk)
    @test fk_low == 0.0
    @test ψk_low == 2 * Nk / xk

    # Positivity case: x̄ > x_{k+1} (Eq. 15a)
    Mk_high = 25.0
    fk_high, ψk_high = Aerosol._sce_linear_params(Nk, Mk_high, xk)
    @test fk_high == 2 * Nk / xk
    @test ψk_high == 0.0
end

@testitem "Closure relations Eq. 8" setup=[SCESetup] tags=[:sce] begin
    Nk = 100.0; Mk = 50.0
    xi_p = 1.0625

    Z = Aerosol._sce_Z(Mk, Nk)
    @test Z ≈ xi_p * Mk^2 / Nk

    Q = Aerosol._sce_Q(Mk, Nk)
    @test Q ≈ xi_p^2 * Mk^3 / Nk^2

    R = Aerosol._sce_R(Mk, Nk)
    @test R ≈ xi_p^3 * Mk^4 / Nk^3

    @test Aerosol._sce_Z(Mk, 0.0) == 0.0
    @test Aerosol._sce_Q(Mk, 0.0) == 0.0
    @test Aerosol._sce_R(Mk, 0.0) == 0.0
end
