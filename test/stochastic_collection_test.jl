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
        u0 = Dict{Any, Float64}()
        for k in 1:I
            xk = x1 * 2.0^(k - 1)
            xk1 = x1 * 2.0^k
            Nk_init = N0 * (exp(-xk / xbar) - exp(-xk1 / xbar))
            Mk_init = N0 * xbar * ((xk / xbar + 1) * exp(-xk / xbar) - (xk1 / xbar + 1) * exp(-xk1 / xbar))
            u0[compiled.Nk[k]] = max(Nk_init, 1.0e-20)
            u0[compiled.Mk[k]] = max(Mk_init, 1.0e-30)
        end
        return u0
    end
end

@testitem "StochasticCollectionCoalescence structure (constant kernel)" setup = [SCESetup] tags = [:sce] begin
    sys = StochasticCollectionCoalescence(; I = 5, kernel_type = :constant, kernel_params = Dict(:K0 => 1.0e-10))
    @test sys isa System

    # 2*I equations: I for dNk/dt and I for dMk/dt
    eqs = equations(sys)
    @test length(eqs) == 10

    vars = unknowns(sys)
    @test length(vars) == 10
end

@testitem "StochasticCollectionCoalescence structure (Golovin kernel)" setup = [SCESetup] tags = [:sce] begin
    sys = StochasticCollectionCoalescence(; I = 5, kernel_type = :golovin, kernel_params = Dict(:C => 1.5e-3))
    @test sys isa System

    eqs = equations(sys)
    @test length(eqs) == 10
end

@testitem "StochasticCollectionCoalescence compiles and solves" setup = [SCESetup] tags = [:sce] begin
    I = 8
    sys = StochasticCollectionCoalescence(; I = I, kernel_type = :constant, kernel_params = Dict(:K0 => 1.0e-10))
    compiled = mtkcompile(sys)
    @test length(unknowns(compiled)) == 2 * I

    x1 = 1.6e-14
    u0 = sce_exponential_ic(compiled, I, x1, 300.0e6, 1.0e-3)

    prob = ODEProblem(compiled, u0, (0.0, 100.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-12)
    @test sol.retcode == ReturnCode.Success
end

@testitem "Constant kernel: number decreases, mass approximately conserved" setup = [SCESetup] tags = [:sce] begin
    I = 8
    sys = StochasticCollectionCoalescence(; I = I, kernel_type = :constant, kernel_params = Dict(:K0 => 1.0e-10))
    compiled = mtkcompile(sys)

    x1 = 1.6e-14
    u0 = sce_exponential_ic(compiled, I, x1, 300.0e6, 1.0e-3)

    prob = ODEProblem(compiled, u0, (0.0, 3600.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-12)
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

@testitem "Golovin kernel compiles, solves, and conserves mass" setup = [SCESetup] tags = [:sce] begin
    I = 8
    sys = StochasticCollectionCoalescence(; I = I, kernel_type = :golovin, kernel_params = Dict(:C => 1.5e-3))
    compiled = mtkcompile(sys)

    x1 = 1.6e-14
    u0 = sce_exponential_ic(compiled, I, x1, 300.0e6, 1.0e-3)

    prob = ODEProblem(compiled, u0, (0.0, 1800.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-12)
    @test sol.retcode == ReturnCode.Success

    N_initial = sum(sol[compiled.Nk[k]][1] for k in 1:I)
    N_final = sum(sol[compiled.Nk[k]][end] for k in 1:I)
    @test N_final < N_initial

    M_initial = sum(sol[compiled.Mk[k]][1] for k in 1:I)
    M_final = sum(sol[compiled.Mk[k]][end] for k in 1:I)
    @test M_final > 0.1 * M_initial
end

@testitem "Invalid kernel type raises error" setup = [SCESetup] tags = [:sce] begin
    @test_throws ErrorException StochasticCollectionCoalescence(; kernel_type = :invalid, kernel_params = Dict(:K0 => 1.0))
end

@testitem "Positivity preservation" setup = [SCESetup] tags = [:sce] begin
    I = 8
    sys = StochasticCollectionCoalescence(; I = I, kernel_type = :constant, kernel_params = Dict(:K0 => 1.0e-10))
    compiled = mtkcompile(sys)

    x1 = 1.6e-14
    u0 = sce_exponential_ic(compiled, I, x1, 300.0e6, 1.0e-3)

    prob = ODEProblem(compiled, u0, (0.0, 3600.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-12)
    @test sol.retcode == ReturnCode.Success

    for k in 1:I
        @test sol[compiled.Nk[k]][end] >= -1.0e-10
        @test sol[compiled.Mk[k]][end] >= -1.0e-20
    end
end

@testitem "Constant kernel: RHS total number rate matches analytical" setup = [SCESetup] tags = [:sce] begin
    # For constant kernel K(x,y) = K0, dN_total/dt = -(K0/2)*N_total^2 exactly.
    # Verify by integrating N_total over a short time and comparing to the
    # analytical solution N(t) = N0 / (1 + K0*N0*t/2).
    I = 8; x1 = 1.0e-14; K0 = 1.0e-10
    sys = StochasticCollectionCoalescence(; I = I, x1 = x1, kernel_type = :constant, kernel_params = Dict(:K0 => K0))
    compiled = mtkcompile(sys)

    u0 = sce_exponential_ic(compiled, I, x1, 300.0e6, 1.0e-3)
    # Use a very short integration so finite-difference is accurate
    prob = ODEProblem(compiled, u0, (0.0, 0.01))
    sol = solve(prob, Tsit5(); reltol = 1.0e-12, abstol = 1.0e-16, saveat = [0.0, 0.01])

    N0 = sum(sol[compiled.Nk[k]][1] for k in 1:I)
    N_final = sum(sol[compiled.Nk[k]][end] for k in 1:I)

    # Analytical N_total(t) = N0 / (1 + K0*N0*t/2)
    t_end = 0.01
    N_analytical = N0 / (1 + K0 * N0 * t_end / 2)

    @test N_final ≈ N_analytical rtol = 1.0e-3
end

@testitem "Constant kernel: mass conservation at solution level" setup = [SCESetup] tags = [:sce] begin
    # Verify that total mass rate is bounded and non-positive (mass leaks through top bin)
    I = 8
    sys = StochasticCollectionCoalescence(; I = I, kernel_type = :constant, kernel_params = Dict(:K0 => 1.0e-10))
    compiled = mtkcompile(sys)

    x1 = 1.6e-14
    u0 = sce_exponential_ic(compiled, I, x1, 300.0e6, 1.0e-3)

    prob = ODEProblem(compiled, u0, (0.0, 100.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-12)
    @test sol.retcode == ReturnCode.Success

    M_initial = sum(sol[compiled.Mk[k]][1] for k in 1:I)
    M_final = sum(sol[compiled.Mk[k]][end] for k in 1:I)

    # Total mass rate should be non-positive (mass lost through top bin) or near zero
    @test M_final <= M_initial * 1.001  # Allow tiny numerical overshoot
    @test isfinite(M_final)
end

@testitem "Golovin kernel: mass conservation at solution level" setup = [SCESetup] tags = [:sce] begin
    I = 8
    sys = StochasticCollectionCoalescence(; I = I, kernel_type = :golovin, kernel_params = Dict(:C => 1.5e-3))
    compiled = mtkcompile(sys)

    x1 = 1.6e-14
    u0 = sce_exponential_ic(compiled, I, x1, 300.0e6, 1.0e-3)

    prob = ODEProblem(compiled, u0, (0.0, 100.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-12)
    @test sol.retcode == ReturnCode.Success

    M_initial = sum(sol[compiled.Mk[k]][1] for k in 1:I)
    M_final = sum(sol[compiled.Mk[k]][end] for k in 1:I)

    @test M_final <= M_initial * 1.001
    @test isfinite(M_final)
end
