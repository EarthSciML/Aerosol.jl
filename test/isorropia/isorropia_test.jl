@testsnippet IsorropiaSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t
    using DynamicQuantities
    using Aerosol
    using OrdinaryDiffEqRosenbrock
    using OrdinaryDiffEqNonlinearSolve
    using NonlinearSolve
    using SciMLBase
end

# =============================================================================
# Equilibrium Constant Tests
# =============================================================================

@testitem "Equilibrium constants at T₀=298.15K" setup = [IsorropiaSetup] tags = [:isorropia] begin
    # At T₀=298.15K, logK_eq should equal logK⁰ = ln(K⁰)
    # Verify K⁰ values match FORTRAN reference (INIT4 in isocom.f)
    fortran_K0 = Dict(
        :r11 => 1.015e-2,   # XK1:  HSO4(aq) ⇌ H(aq) + SO4(aq)
        :r12 => 5.764e1,    # XK21: NH3(g) ⇌ NH3(aq)
        :r13 => 1.805e-5,   # XK22: NH3(aq) ⇌ NH4(aq) + OH(aq)
        :r14 => 2.511e6,    # XK4:  HNO3(g) ⇌ H(aq) + NO3(aq)
        :r15 => 2.1e5,      # XK41: HNO3(g) ⇌ HNO3(aq)
        :r16 => 1.971e6,    # XK3:  HCl(g) ⇌ H(aq) + Cl(aq)
        :r17 => 2.5e3,      # XK31: HCl(g) ⇌ HCl(aq)
        :r18 => 1.01e-14,   # XKW:  H2O ⇌ H(aq) + OH(aq)
        :r19 => 4.799e-1,   # XK5:  Na2SO4(s) ⇌ 2Na(aq) + SO4(aq)
        :r20 => 1.817e0,    # XK7:  (NH4)2SO4(s) ⇌ 2NH4(aq) + SO4(aq)
        :r21 => 1.086e-16,  # XK6:  NH4Cl(s) ⇌ NH3(g) + HCl(g)
        :r22 => 1.197e1,    # XK9:  NaNO3(s) ⇌ Na(aq) + NO3(aq)
        :r23 => 3.766e1,    # XK8:  NaCl(s) ⇌ Na(aq) + Cl(aq)
        :r24 => 2.413e4,    # XK11: NaHSO4(s) ⇌ Na(aq) + HSO4(aq)
        :r25 => 4.199e-17,  # XK10: NH4NO3(s) ⇌ NH3(g) + HNO3(g)
        :r26 => 1.382e2,    # XK12: NH4HSO4(s) ⇌ NH4(aq) + HSO4(aq)
        :r27 => 2.972e1,    # XK13: (NH4)3H(SO4)2(s) ⇌ 3NH4(aq) + HSO4(aq) + SO4(aq)
        :r4 => 1.569e-2,   # XK17: K2SO4(s) ⇌ 2K(aq) + SO4(aq)
        :r5 => 24.016,     # XK18: KHSO4(s) ⇌ K(aq) + HSO4(aq)
        :r6 => 0.872,      # XK19: KNO3(s) ⇌ K(aq) + NO3(aq)
        :r7 => 8.68,       # XK20: KCl(s) ⇌ K(aq) + Cl(aq)
        :r8 => 1.079e5,    # XK23: MgSO4(s) ⇌ Mg(aq) + SO4(aq)
        :r9 => 2.507e15,   # XK24: Mg(NO3)2(s) ⇌ Mg(aq) + 2NO3(aq)
        :r10 => 9.557e21,   # XK25: MgCl2(s) ⇌ Mg(aq) + 2Cl(aq)
    )

    eq = ISORROPIA.EquilibriumConstants()
    # Wrap EquilibriumConstants to constrain T (a variable) with a parameter
    @parameters T_val = 298.15 [unit = u"K"]
    wrapper = System([eq.T ~ T_val], t; systems = [eq], name = :test_eq)
    sys = mtkcompile(wrapper)
    # At T=298.15, equilibrium constants should equal K⁰
    prob = NonlinearProblem(sys, [])
    sol = solve(prob, NewtonRaphson())
    @test sol.retcode == SciMLBase.ReturnCode.Success

    for (rname, K0) in fortran_K0
        rsym = getproperty(eq, rname)
        logK_val = sol[rsym.logK_eq]
        K_computed = exp(logK_val)
        @test K_computed ≈ K0 rtol = 1.0e-3
    end
end

@testitem "Equilibrium constants temperature dependence" setup = [IsorropiaSetup] tags =
    [:isorropia] begin
    # Verify temperature correction against FORTRAN reference at T=310K
    # FORTRAN formula: K(T) = K0 * exp(A*(T0/T - 1) + B*(1 + ln(T0/T) - T0/T))
    T = 310.0
    T0 = 298.15
    T0T = T0 / T
    COEF = 1.0 + log(T0T) - T0T

    # FORTRAN reference values: (K0, A, B)
    fortran_params = Dict(
        :r11 => (1.015e-2, 8.85, 25.14),       # XK1
        :r12 => (5.764e1, 13.79, -5.393),       # XK21
        :r13 => (1.805e-5, -1.5, 26.92),       # XK22
        :r14 => (2.511e6, 29.17, 16.83),        # XK4
        :r16 => (1.971e6, 30.2, 19.91),        # XK3
        :r18 => (1.01e-14, -22.52, 26.92),      # XKW
        :r19 => (4.799e-1, 0.98, 39.5),        # XK5 (note: FORTRAN has 39.5, paper may have 39.75)
        :r20 => (1.817e0, -2.65, 38.57),        # XK7
        :r21 => (1.086e-16, -71.0, 2.4),       # XK6
        :r22 => (1.197e1, -8.22, 16.01),        # XK9
        :r23 => (3.766e1, -1.56, 16.9),        # XK8
        :r24 => (2.413e4, 0.79, 14.746),        # XK11
        :r26 => (1.382e2, -2.87, 15.83),        # XK12
        :r27 => (2.972e1, -5.19, 54.4),        # XK13
        :r4 => (1.569e-2, -9.585, 45.81),      # XK17
        :r5 => (24.016, -8.423, 17.96),        # XK18
        :r6 => (0.872, -14.08, 19.39),         # XK19 (note: FORTRAN uses -14.08, Julia uses -14.075)
        :r7 => (8.68, -6.902, 19.95),          # XK20 (note: value diffs from Julia's -6.167)
    )

    eq = ISORROPIA.EquilibriumConstants()
    # Wrap EquilibriumConstants to constrain T (a variable) with a parameter
    @parameters T_val = T [unit = u"K"]
    wrapper = System([eq.T ~ T_val], t; systems = [eq], name = :test_eq)
    sys = mtkcompile(wrapper)
    prob = NonlinearProblem(sys, [])
    sol = solve(prob, NewtonRaphson())
    @test sol.retcode == SciMLBase.ReturnCode.Success

    for (rname, (K0, A, B)) in fortran_params
        # Compute expected K using FORTRAN formula
        K_expected = K0 * exp(A * (T0T - 1.0) + B * COEF)

        rsym = getproperty(eq, rname)
        K_computed = exp(sol[rsym.logK_eq])

        # Allow some tolerance for reactions where Julia's H/C values differ slightly from FORTRAN
        @test K_computed ≈ K_expected rtol = 0.05
    end
end

# =============================================================================
# Model Compilation Tests
# =============================================================================

@testitem "Isorropia model compilation" setup = [IsorropiaSetup] tags = [:isorropia] begin
    isrpa = Isorropia()
    sys = mtkcompile(isrpa)
    @test sys isa System

    # Verify key variables exist
    @test hasproperty(sys, :aq)
    @test hasproperty(sys, :g)
    @test hasproperty(sys, :eq)
    @test hasproperty(sys, :NH)
    @test hasproperty(sys, :Na)
    @test hasproperty(sys, :Ca)
    @test hasproperty(sys, :K)
    @test hasproperty(sys, :Mg)
    @test hasproperty(sys, :Cl)
    @test hasproperty(sys, :NO3)
    @test hasproperty(sys, :SO4)
end

@testitem "Isorropia aerosol type classification" setup = [IsorropiaSetup] tags =
    [:isorropia] begin
    # Test the aerosol type smoothed classification based on sulfate ratios (Section 3.1)
    # Type 1: R₁ < 1 (sulfate rich, free acid)
    # Type 2: 1 ≤ R₁ < 2 (sulfate rich)
    # Type 3: R₁ ≥ 2, R₂ < 2 (sulfate poor, crustal & sodium poor)
    # Type 4: R₁ ≥ 2, R₂ ≥ 2, R₃ < 2 (sulfate poor, crustal & sodium rich, crustal poor)
    # Type 5: R₁ ≥ 2, R₂ ≥ 2, R₃ ≥ 2 (sulfate poor, crustal & sodium rich, crustal rich)

    # Use the smooth tanh approximations directly
    function type_fractions(R1, R2, R3)
        s(x) = (tanh(x * 30) + 1) / 2
        t1 = 1 - s(R1 - 1)
        t2 = min(s(R1 - 1), 1 - s(R1 - 2))
        t3 = min(s(R1 - 2), 1 - s(R2 - 2))
        t4 = min(s(R1 - 2), s(R2 - 2), 1 - s(R3 - 2))
        t5 = min(s(R1 - 2), s(R2 - 2), s(R3 - 2))
        return (t1, t2, t3, t4, t5)
    end

    # Type 1 dominant: R₁ = 0.5, R₂ = 0.3, R₃ = 0.1
    tf = type_fractions(0.5, 0.3, 0.1)
    @test tf[1] ≈ 1.0 atol = 0.01
    @test sum(tf) ≈ 1.0 atol = 0.01

    # Type 3 dominant: R₁ = 5, R₂ = 0.5, R₃ = 0.2
    tf = type_fractions(5.0, 0.5, 0.2)
    @test tf[3] ≈ 1.0 atol = 0.01
    @test sum(tf) ≈ 1.0 atol = 0.01

    # Type 4 dominant: R₁ = 5, R₂ = 5, R₃ = 0.5
    tf = type_fractions(5.0, 5.0, 0.5)
    @test tf[4] ≈ 1.0 atol = 0.01
    @test sum(tf) ≈ 1.0 atol = 0.01

    # Type 5 dominant: R₁ = 5, R₂ = 5, R₃ = 5
    tf = type_fractions(5.0, 5.0, 5.0)
    @test tf[5] ≈ 1.0 atol = 0.01
    @test sum(tf) ≈ 1.0 atol = 0.01

    # Verify type4 ≠ type5 (this was a bug: they were identical)
    tf_mid = type_fractions(5.0, 5.0, 2.0)
    @test tf_mid[4] ≈ tf_mid[5] atol = 0.1  # At exactly R₃=2, both should be ~0.5
    tf_low_r3 = type_fractions(5.0, 5.0, 0.5)
    tf_high_r3 = type_fractions(5.0, 5.0, 3.5)
    @test tf_low_r3[4] > 0.9  # Low R₃ → type 4
    @test tf_low_r3[5] < 0.1
    @test tf_high_r3[5] > 0.9  # High R₃ → type 5
    @test tf_high_r3[4] < 0.1
end

@testitem "Isorropia ODE solve" setup = [IsorropiaSetup] tags = [:isorropia] begin
    isrpa = Isorropia()
    sys = mtkcompile(isrpa)

    prob = ODEProblem(
        sys,
        [sys.RH => 0.7],
        (0.0, 10.0);
        initializealg = BrownFullBasicInit(nlsolve = RobustMultiNewton()),
    )

    sol = solve(prob, Rosenbrock23())

    # Check that totals are conserved (D(total) ~ 0 for all species)
    for sp in [
            sys.NH.total,
            sys.Na.total,
            sys.Ca.total,
            sys.K.total,
            sys.Mg.total,
            sys.Cl.total,
            sys.NO3.total,
            sys.SO4.total,
        ]
        vals = sol[sp]
        @test vals[1] ≈ vals[end] rtol = 1.0e-6
    end
end

@testitem "Salt group helper functions" setup = [IsorropiaSetup] tags = [:isorropia] begin
    isrpa = Isorropia()

    # Test salt_group returns correct number of salts for each group
    @test length(ISORROPIA.salt_group(isrpa.aq, :NH4, :M)) == 5
    @test length(ISORROPIA.salt_group(isrpa.aq, :Na, :M)) == 4
    @test length(ISORROPIA.salt_group(isrpa.aq, :Ca, :M)) == 3
    @test length(ISORROPIA.salt_group(isrpa.aq, :K, :M)) == 4
    @test length(ISORROPIA.salt_group(isrpa.aq, :Mg, :M)) == 3
    @test length(ISORROPIA.salt_group(isrpa.aq, :Cl, :M)) == 5
    @test length(ISORROPIA.salt_group(isrpa.aq, :NO3, :M)) == 5
    @test length(ISORROPIA.salt_group(isrpa.aq, :SO4, :M)) == 5

    # Test salt_group_ν returns correct type
    @test ISORROPIA.salt_group_ν(:NH4) == :ν_cation
    @test ISORROPIA.salt_group_ν(:Cl) == :ν_anion
    @test ISORROPIA.salt_group_ν(:SO4) == :ν_anion
end
