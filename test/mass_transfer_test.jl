@testsnippet MassTransferSetup begin
    using Test
    using ModelingToolkit
    using NonlinearSolve
    using Aerosol
end

@testitem "MeanMolecularSpeed structure" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = MeanMolecularSpeed()
    @test sys isa System
    compiled = mtkcompile(sys)
    # Algebraic system with no dynamics has no unknowns after compilation
    @test length(equations(sys)) == 1
end

@testitem "MeanMolecularSpeed values" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    # Test mean molecular speed for air at 298 K
    # c̄ = √(8kT/(πm)) ≈ 464 m/s for air at 298 K
    sys = MeanMolecularSpeed()
    compiled = mtkcompile(sys)

    # Set parameters for air (M_A = 0.029 kg/mol)
    prob = NonlinearProblem(compiled, [], [compiled.T => 298.15, compiled.M_A => 0.029])
    sol = solve(prob)

    c_bar = sol[compiled.c_bar]
    # Expected value for air at 298 K: approximately 464 m/s
    @test c_bar ≈ 464.0 rtol=0.05
end

@testitem "KnudsenNumber structure" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = KnudsenNumber()
    @test sys isa System
    @test length(equations(sys)) == 1
end

@testitem "KnudsenNumber values" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = KnudsenNumber()
    compiled = mtkcompile(sys)

    # Test with λ = 65 nm, R_p = 100 nm → Kn = 0.65
    prob = NonlinearProblem(compiled, [],
        [compiled.λ => 6.5e-8, compiled.R_p => 1.0e-7])
    sol = solve(prob)

    Kn = sol[compiled.Kn]
    @test Kn ≈ 0.65 rtol=1e-6
end

@testitem "FuchsSutugin structure" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = FuchsSutugin()
    @test sys isa System
    @test length(equations(sys)) == 1
end

@testitem "FuchsSutugin limiting behavior" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = FuchsSutugin()
    compiled = mtkcompile(sys)

    # Continuum limit: Kn → 0, f_FS → 1 (for α = 1)
    prob_cont = NonlinearProblem(compiled, [],
        [compiled.Kn => 0.001, compiled.α => 1.0])
    sol_cont = solve(prob_cont)
    @test sol_cont[compiled.f_FS] ≈ 1.0 rtol=0.01

    # Kinetic limit: Kn → ∞, f_FS → 0.75α·Kn⁻¹ (approaches 0)
    prob_kin = NonlinearProblem(compiled, [],
        [compiled.Kn => 100.0, compiled.α => 1.0])
    sol_kin = solve(prob_kin)
    @test sol_kin[compiled.f_FS] < 0.01

    # Transition regime: Kn = 1
    prob_trans = NonlinearProblem(compiled, [],
        [compiled.Kn => 1.0, compiled.α => 1.0])
    sol_trans = solve(prob_trans)
    f_trans = sol_trans[compiled.f_FS]
    # At Kn = 1, α = 1: f_FS ≈ 0.75(2)/(1+1+0.283+0.75) ≈ 0.494
    @test f_trans ≈ 0.494 rtol=0.01
end

@testitem "FuchsSutugin accommodation effect" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    # Test that lower accommodation coefficient reduces mass transfer
    sys = FuchsSutugin()
    compiled = mtkcompile(sys)

    Kn_val = 1.0

    prob_α1 = NonlinearProblem(compiled, [],
        [compiled.Kn => Kn_val, compiled.α => 1.0])
    sol_α1 = solve(prob_α1)

    prob_α01 = NonlinearProblem(compiled, [],
        [compiled.Kn => Kn_val, compiled.α => 0.1])
    sol_α01 = solve(prob_α01)

    # Lower α should give lower f_FS
    @test sol_α01[compiled.f_FS] < sol_α1[compiled.f_FS]
end

@testitem "Dahneke structure" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = Dahneke()
    @test sys isa System
    @test length(equations(sys)) == 1
end

@testitem "Dahneke limiting behavior" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = Dahneke()
    compiled = mtkcompile(sys)

    # Continuum limit: Kn → 0, f_D → 1 (for α = 1)
    prob_cont = NonlinearProblem(compiled, [],
        [compiled.Kn => 0.001, compiled.α => 1.0])
    sol_cont = solve(prob_cont)
    @test sol_cont[compiled.f_D] ≈ 1.0 rtol=0.01
end

@testitem "MaxwellianFlux structure" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = MaxwellianFlux()
    @test sys isa System
    @test length(equations(sys)) == 1
end

@testitem "MaxwellianFlux values" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = MaxwellianFlux()
    compiled = mtkcompile(sys)

    # Test J_c = 4πR_p D_g (c_∞ - c_s)
    R_p = 1.0e-6  # 1 μm
    D_g = 2.0e-5  # m²/s
    c_inf = 1.0e-6  # mol/m³
    c_s = 0.0

    expected_J_c = 4 * π * R_p * D_g * (c_inf - c_s)

    prob = NonlinearProblem(compiled, [],
        [compiled.R_p => R_p, compiled.D_g => D_g,
         compiled.c_inf => c_inf, compiled.c_s => c_s])
    sol = solve(prob)

    @test sol[compiled.J_c] ≈ expected_J_c rtol=1e-10
end

@testitem "MassTransferCoefficient structure" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = MassTransferCoefficient()
    @test sys isa System
end

@testitem "MassTransfer comprehensive" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = MassTransfer()
    @test sys isa System
    compiled = mtkcompile(sys)

    # Check that all expected equations exist in the original system
    eq_strs = string.(equations(sys))
    @test any(occursin("c_bar", s) for s in eq_strs)
    @test any(occursin("Kn", s) for s in eq_strs)
    @test any(occursin("f_FS", s) for s in eq_strs)
    @test any(occursin("J_c", s) for s in eq_strs)
end

@testitem "MassTransfer regime dependence" setup=[MassTransferSetup] tags=[:mass_transfer] begin
    sys = MassTransfer()
    compiled = mtkcompile(sys)

    # Large particle (continuum regime, Kn ≪ 1): J ≈ J_c
    prob_large = NonlinearProblem(compiled, [],
        [compiled.R_p => 1.0e-5, compiled.α => 1.0,
         compiled.c_inf => 1.0e-6, compiled.c_s => 0.0])
    sol_large = solve(prob_large)
    J_large = sol_large[compiled.J]
    J_c_large = sol_large[compiled.J_c]
    # In continuum regime, J should be close to J_c
    @test J_large / J_c_large > 0.9

    # Small particle (kinetic regime, Kn ≫ 1): J < J_c
    prob_small = NonlinearProblem(compiled, [],
        [compiled.R_p => 1.0e-8, compiled.α => 1.0,
         compiled.c_inf => 1.0e-6, compiled.c_s => 0.0])
    sol_small = solve(prob_small)
    J_small = sol_small[compiled.J]
    J_c_small = sol_small[compiled.J_c]
    # In kinetic regime, J should be noticeably less than J_c
    # At R_p = 10 nm, Kn is about 6-7, so J/J_c ≈ 0.1
    @test J_small / J_c_small < 0.15
end
