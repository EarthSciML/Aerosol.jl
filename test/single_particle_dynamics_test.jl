@testsnippet SPDSetup begin
    using Test
    using ModelingToolkit
    using Aerosol
end

@testitem "MeanFreePath structural" setup=[SPDSetup] tags=[:spd] begin
    sys=MeanFreePath()
    # Check that we have the equation
    @test length(equations(sys)) >= 1
    # Check variables exist
    vars=unknowns(sys)
    @test any(v -> contains(string(v), "λ"), vars)
end

@testitem "MeanFreePath at standard conditions" setup=[SPDSetup] tags=[:spd] begin
    # Eq. 9.7: λ_air(298K, 1atm) = 0.0651 μm = 6.51e-8 m
    sys=MeanFreePath()
    csys=mtkcompile(sys)

    # Default parameters are T=298.15K, P=101325Pa, μ=1.8e-5 kg/(m·s)
    # For purely algebraic systems, we can get the value from defaults
    prob=NonlinearProblem(csys, Dict(), Dict())
    sol=solve(prob)

    λ_computed=sol[csys.λ]
    λ_expected=6.51e-8  # m (from Table 9.7)

    @test λ_computed ≈ λ_expected rtol=0.05
end

@testitem "SlipCorrection structural" setup=[SPDSetup] tags=[:spd] begin
    sys=SlipCorrection()
    @test length(equations(sys)) >= 2  # C_c and Kn equations
    vars=unknowns(sys)
    @test any(v -> contains(string(v), "C_c"), vars)
    @test any(v -> contains(string(v), "Kn"), vars)
end

@testitem "SlipCorrection values from Table 9.3" setup=[SPDSetup] tags=[:spd] begin
    # Table 9.3 validation: Slip correction factors at 298K, 1atm
    # λ = 0.0651 μm = 6.51e-8 m

    test_cases=[
        # (D_p in m, expected C_c)
        (1e-8, 22.2),   # 0.01 μm
        (1e-7, 2.85),   # 0.1 μm
        (1e-6, 1.164),  # 1 μm
        (1e-5, 1.016)  # 10 μm
    ]

    λ_air=6.51e-8  # m

    for (D_p, C_c_expected) in test_cases
        sys=SlipCorrection()
        csys=mtkcompile(sys)

        prob=NonlinearProblem(csys, Dict(), Dict(csys.D_p=>D_p, csys.λ=>λ_air))
        sol=solve(prob)

        C_c_computed=sol[csys.C_c]

        @test C_c_computed ≈ C_c_expected rtol=0.05
    end
end

@testitem "SlipCorrection Knudsen number" setup=[SPDSetup] tags=[:spd] begin
    sys=SlipCorrection()
    csys=mtkcompile(sys)

    D_p=1e-7  # 0.1 μm
    λ=6.51e-8  # m

    prob=NonlinearProblem(csys, Dict(), Dict(csys.D_p=>D_p, csys.λ=>λ))
    sol=solve(prob)

    # Kn = 2λ/D_p = 2 * 6.51e-8 / 1e-7 = 1.302
    Kn_expected=2*λ/D_p
    @test sol[csys.Kn] ≈ Kn_expected rtol=1e-6
end

@testitem "SettlingVelocity structural" setup=[SPDSetup] tags=[:spd] begin
    sys=SettlingVelocity()
    @test length(equations(sys)) >= 2  # v_t and τ equations
end

@testitem "SettlingVelocity relaxation time from Table 9.4" setup=[SPDSetup] tags=[:spd] begin
    # Table 9.4: Relaxation times for ρ_p = 1 g/cm³ = 1000 kg/m³
    # at 298K, 1atm

    test_cases=[
        # (D_p in m, C_c, τ_expected in s)
        (1e-6, 1.164, 3.6e-6),    # 1 μm
        (1e-5, 1.016, 3.14e-4)   # 10 μm
    ]

    for (D_p, C_c, τ_expected) in test_cases
        sys=SettlingVelocity()
        csys=mtkcompile(sys)

        prob=NonlinearProblem(csys, Dict(), Dict(
            csys.D_p=>D_p,
            csys.ρ_p=>1000.0,
            csys.C_c=>C_c
        ))
        sol=solve(prob)

        @test sol[csys.τ] ≈ τ_expected rtol=0.1
    end
end

@testitem "SettlingVelocity v_t = τg" setup=[SPDSetup] tags=[:spd] begin
    # Terminal velocity should equal τ * g
    sys=SettlingVelocity()
    csys=mtkcompile(sys)

    D_p=1e-6
    ρ_p=1000.0
    C_c=1.164

    prob=NonlinearProblem(csys, Dict(), Dict(
        csys.D_p=>D_p,
        csys.ρ_p=>ρ_p,
        csys.C_c=>C_c
    ))
    sol=solve(prob)

    g=9.807
    @test sol[csys.v_t] ≈ sol[csys.τ] * g rtol=1e-6
end

@testitem "BrownianDiffusion structural" setup=[SPDSetup] tags=[:spd] begin
    sys=BrownianDiffusion()
    @test length(equations(sys)) >= 5  # D_B, B, c_p, λ_p, m_p equations
end

@testitem "BrownianDiffusion coefficient from Table 9.5" setup=[SPDSetup] tags=[:spd] begin
    # Table 9.5: Brownian motion quantities at 298K, 1atm

    test_cases=[
        # (D_p in m, C_c, D_expected in m²/s)
        (1e-8, 22.2, 5.38e-8),    # 0.01 μm: D = 5.38e-4 cm²/s = 5.38e-8 m²/s
        (1e-7, 2.85, 6.91e-10),   # 0.1 μm: D = 6.91e-6 cm²/s = 6.91e-10 m²/s
        (1e-6, 1.164, 2.82e-11)  # 1 μm: D = 2.82e-7 cm²/s = 2.82e-11 m²/s
    ]

    for (D_p, C_c, D_expected) in test_cases
        sys=BrownianDiffusion()
        csys=mtkcompile(sys)

        prob=NonlinearProblem(csys, Dict(), Dict(
            csys.D_p=>D_p,
            csys.ρ_p=>1000.0,
            csys.T=>298.0,
            csys.C_c=>C_c
        ))
        sol=solve(prob)

        @test sol[csys.D_B] ≈ D_expected rtol=0.15
    end
end

@testitem "BrownianDiffusion Einstein relation D = BkT" setup=[SPDSetup] tags=[:spd] begin
    # Eq. 9.79: D = B k T
    sys=BrownianDiffusion()
    csys=mtkcompile(sys)

    T=298.0
    k_B=1.381e-23

    prob=NonlinearProblem(csys, Dict(), Dict(csys.T=>T))
    sol=solve(prob)

    D_from_BkT=sol[csys.B]*k_B*T
    @test sol[csys.D_B] ≈ D_from_BkT rtol=1e-6
end

@testitem "ElectricalMobility structural" setup=[SPDSetup] tags=[:spd] begin
    sys=ElectricalMobility()
    @test length(equations(sys)) >= 3  # q, B_e, v_e equations
end

@testitem "ElectricalMobility v_e = B_e * E" setup=[SPDSetup] tags=[:spd] begin
    # Eq. 9.51: v_e = B_e E
    sys=ElectricalMobility()
    csys=mtkcompile(sys)

    E=5000.0  # V/m

    prob=NonlinearProblem(csys, Dict(), Dict(csys.E=>E))
    sol=solve(prob)

    @test sol[csys.v_e] ≈ sol[csys.B_e] * E rtol=1e-6
end

@testitem "StokesNumber structural" setup=[SPDSetup] tags=[:spd] begin
    sys=StokesNumber()
    @test length(equations(sys)) >= 3  # τ, s_p, St equations
end

@testitem "StokesNumber s_p = τU" setup=[SPDSetup] tags=[:spd] begin
    # Eq. 9.98: s_p = τ U
    sys=StokesNumber()
    csys=mtkcompile(sys)

    U=10.0  # m/s

    prob=NonlinearProblem(csys, Dict(), Dict(csys.U=>U))
    sol=solve(prob)

    @test sol[csys.s_p] ≈ sol[csys.τ] * U rtol=1e-6
end

@testitem "StokesNumber St = τU/L" setup=[SPDSetup] tags=[:spd] begin
    # Eq. 9.101: St = τ u_0 / L
    sys=StokesNumber()
    csys=mtkcompile(sys)

    U=10.0  # m/s
    L=0.001  # m

    prob=NonlinearProblem(csys, Dict(), Dict(csys.U=>U, csys.L=>L))
    sol=solve(prob)

    @test sol[csys.St] ≈ sol[csys.τ] * U / L rtol=1e-6
end

@testitem "AerodynamicDiameter structural" setup=[SPDSetup] tags=[:spd] begin
    sys=AerodynamicDiameter()
    @test length(equations(sys)) >= 2  # D_ca, D_ca_continuum equations
end

@testitem "AerodynamicDiameter continuum limit" setup=[SPDSetup] tags=[:spd] begin
    # Eq. 9.112: D_ca/D_p = sqrt(ρ_p/ρ_p°) in continuum limit
    # Example from text: ρ_p = 1.5 g/cm³, D_p = 1 μm -> D_ca = 1.242 μm (with slip)
    # In continuum: D_ca = D_p * sqrt(1.5) = 1.225 μm

    sys=AerodynamicDiameter()
    csys=mtkcompile(sys)

    D_p=1e-6  # 1 μm
    ρ_p=1500.0  # 1.5 g/cm³

    prob=NonlinearProblem(csys, Dict(), Dict(csys.D_p=>D_p, csys.ρ_p=>ρ_p))
    sol=solve(prob)

    D_ca_expected=D_p*sqrt(ρ_p/1000.0)
    @test sol[csys.D_ca_continuum] ≈ D_ca_expected rtol=1e-6
    @test sol[csys.D_ca_continuum] ≈ 1.225e-6 rtol=0.01
end

@testitem "SingleParticleDynamics structural" setup=[SPDSetup] tags=[:spd] begin
    sys=SingleParticleDynamics()
    # Check we have all the expected equations (11 variables)
    @test length(equations(sys)) >= 11
end

@testitem "SingleParticleDynamics comprehensive validation" setup=[SPDSetup] tags=[:spd] begin
    # Test the comprehensive component at standard conditions
    sys=SingleParticleDynamics()
    csys=mtkcompile(sys)

    # 1 μm particle at standard conditions
    D_p=1e-6
    ρ_p=1000.0

    prob=NonlinearProblem(csys, Dict(), Dict(csys.D_p=>D_p, csys.ρ_p=>ρ_p))
    sol=solve(prob)

    # Validate against Table 9.3: C_c ≈ 1.164 for 1 μm
    @test sol[csys.C_c] ≈ 1.164 rtol=0.05

    # Validate mean free path (Eq. 9.7)
    @test sol[csys.λ] ≈ 6.51e-8 rtol=0.05

    # Validate Knudsen number
    Kn_expected=2*sol[csys.λ]/D_p
    @test sol[csys.Kn] ≈ Kn_expected rtol=1e-6

    # Validate relaxation time (Table 9.4): τ ≈ 3.6e-6 s for 1 μm
    @test sol[csys.τ] ≈ 3.6e-6 rtol=0.1

    # Validate v_t = τ * g
    @test sol[csys.v_t] ≈ sol[csys.τ] * 9.807 rtol=1e-6

    # Validate diffusion coefficient (Table 9.5): D ≈ 2.82e-11 m²/s for 1 μm
    @test sol[csys.D_B] ≈ 2.82e-11 rtol=0.15
end

@testitem "SingleParticleDynamics temperature dependence" setup=[SPDSetup] tags=[:spd] begin
    # Test that mean free path increases with temperature (via viscosity effect)
    sys=SingleParticleDynamics()
    csys=mtkcompile(sys)

    # At 298K
    prob1=NonlinearProblem(csys, Dict(), Dict(csys.T=>298.0))
    sol1=solve(prob1)

    # At higher temperature (320K)
    prob2=NonlinearProblem(csys, Dict(), Dict(csys.T=>320.0))
    sol2=solve(prob2)

    # Mean free path should increase with temperature
    @test sol2[csys.λ] > sol1[csys.λ]

    # Diffusion coefficient should increase with temperature
    @test sol2[csys.D_B] > sol1[csys.D_B]
end

@testitem "SingleParticleDynamics size dependence" setup=[SPDSetup] tags=[:spd] begin
    # Test qualitative size dependence
    sys=SingleParticleDynamics()
    csys=mtkcompile(sys)

    # Small particle (0.1 μm)
    prob1=NonlinearProblem(csys, Dict(), Dict(csys.D_p=>1e-7))
    sol1=solve(prob1)

    # Large particle (10 μm)
    prob2=NonlinearProblem(csys, Dict(), Dict(csys.D_p=>1e-5))
    sol2=solve(prob2)

    # Smaller particles have larger slip correction
    @test sol1[csys.C_c] > sol2[csys.C_c]

    # Smaller particles have larger diffusion coefficient
    @test sol1[csys.D_B] > sol2[csys.D_B]

    # Larger particles have larger settling velocity
    @test sol2[csys.v_t] > sol1[csys.v_t]

    # Larger particles have larger relaxation time
    @test sol2[csys.τ] > sol1[csys.τ]
end
