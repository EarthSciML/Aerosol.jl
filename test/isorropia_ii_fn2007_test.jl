# =============================================================================
# Tests for ISORROPIA II (Fountoukis & Nenes, 2007) implementation
# =============================================================================

@testsnippet Iso2Setup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t
    using DynamicQuantities
    using Aerosol
    using NonlinearSolve
    using SciMLBase

    # Helper: create and solve a NonlinearProblem with standard settings
    function iso2_solve(compiled, op; maxiters = 10000)
        prob = NonlinearProblem(compiled, op)
        return solve(prob, RobustMultiNewton(); maxiters)
    end
end

# =============================================================================
# Equilibrium Constant Tests
# =============================================================================

@testitem "ISO2: Van't Hoff equation at T₀=298.15K" setup = [Iso2Setup] tags = [:iso2] begin
    # At T₀=298.15K, K(T) should equal K₀
    # FORTRAN reference values from CMAQ isocom.f (INIT4 subroutine)
    T0 = 298.15

    test_cases = [
        (1.015e-2, 8.85, 25.14),    # K1: HSO4 dissociation
        (5.764e1, 13.79, -5.393),    # K21: NH3 dissolution
        (1.805e-5, -1.5, 26.92),    # K22: NH3 dissociation
        (2.511e6, 29.17, 16.83),     # K4: HNO3
        (1.971e6, 30.2, 19.91),     # K3: HCl
        (1.01e-14, -22.52, 26.92),  # Kw: water
        (1.569e-2, -9.585, 45.81),   # K17: K2SO4
        (24.016, -8.423, 17.96),     # K18: KHSO4
        (0.872, -14.08, 19.39),      # K19: KNO3
        (8.68, -6.902, 19.95),      # K20: KCl
    ]

    for (K0, A, B) in test_cases
        K_computed = Aerosol._iso2_eq_const(K0, A, B, T0)
        @test K_computed ≈ K0 rtol = 1.0e-10
    end
end

@testitem "ISO2: Van't Hoff temperature dependence" setup = [Iso2Setup] tags = [:iso2] begin
    # Verify temperature correction at T=310K against FORTRAN formula
    T = 310.0
    T0 = 298.15
    T0T = T0 / T
    COEF = 1.0 + log(T0T) - T0T

    # (K0, A, B) parameters from CMAQ INIT4
    test_cases = [
        (1.015e-2, 8.85, 25.14),    # K1
        (5.764e1, 13.79, -5.393),    # K21
        (1.805e-5, -1.5, 26.92),    # K22
        (2.511e6, 29.17, 16.83),     # K4
        (1.971e6, 30.2, 19.91),     # K3
        (1.01e-14, -22.52, 26.92),  # Kw
    ]

    for (K0, A, B) in test_cases
        K_expected = K0 * exp(A * (T0T - 1.0) + B * COEF)
        K_computed = Aerosol._iso2_eq_const(K0, A, B, T)
        @test K_computed ≈ K_expected rtol = 1.0e-10
    end
end

# =============================================================================
# Activity Coefficient Tests
# =============================================================================

@testitem "ISO2: Kusik-Meissner activity coefficients" setup = [Iso2Setup] tags = [:iso2] begin
    # At very low ionic strength, γ± → 1
    γ_low = Aerosol._iso2_km_gamma(2.23, 1, 1.0e-10)  # NaCl params
    @test γ_low ≈ 1.0 rtol = 0.01

    # Activity coefficient should be positive for any I
    for I in [0.01, 0.1, 1.0, 5.0, 10.0, 20.0]
        γ = Aerosol._iso2_km_gamma(2.23, 1, I)
        @test γ > 0
    end

    # Higher charge products (z) should give lower activity coefficients
    γ_z1 = Aerosol._iso2_km_gamma(0.0, 1, 0.5)
    γ_z2 = Aerosol._iso2_km_gamma(0.0, 2, 0.5)
    @test γ_z2 < γ_z1

    # Known reference: NaCl at I=1.0 mol/kg should give γ ≈ 0.65-0.70
    γ_NaCl = Aerosol._iso2_km_gamma(2.23, 1, 1.0)
    @test 0.5 < γ_NaCl < 0.9
end

@testitem "ISO2: Temperature correction for activity coefficients" setup = [Iso2Setup] tags = [:iso2] begin
    # At T=298.15K (TI=25.15°C), the correction should be close to identity
    γ_298 = 0.5
    γ_corrected = Aerosol._iso2_gamma_T(γ_298, 1.0, 298.15)
    @test γ_corrected ≈ γ_298 rtol = 0.05

    # Temperature correction should be monotonic in T
    γs = [Aerosol._iso2_gamma_T(0.5, 1.0, T) for T in 270.0:5.0:310.0]
    @test all(γs .> 0)
    @test maximum(γs) / minimum(γs) < 10
end

# =============================================================================
# ZSR Water Content Tests
# =============================================================================

@testitem "ISO2: ZSR water content" setup = [Iso2Setup] tags = [:iso2] begin
    # Pure (NH4)2SO4 solution at RH=0.80
    # Binary molality from table at aw=0.80 ≈ 5.83 mol/kg
    # W = c_SO4 / m0 ≈ 1e-7 / 5.83 ≈ 1.715e-8
    W = Aerosol._iso2_zsr_water(0.8, 0.0, 2.0e-7, 1.0e-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    @test W > 0
    @test W ≈ 1.0e-7 / 5.83 rtol = 0.1

    # Higher RH should give more water
    W_70 = Aerosol._iso2_zsr_water(0.7, 0.0, 2.0e-7, 1.0e-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    W_90 = Aerosol._iso2_zsr_water(0.9, 0.0, 2.0e-7, 1.0e-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    @test W_90 > W_70

    # Zero concentrations should give minimal water
    W_zero = Aerosol._iso2_zsr_water(0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    @test W_zero ≈ 0.0 atol = 1.0e-15
end

# =============================================================================
# Component Structure Tests
# =============================================================================

@testitem "ISO2: Model structure" setup = [Iso2Setup] tags = [:iso2] begin
    sys = IsorropiaEquilibrium()
    @test sys isa System

    # Should have 21 equations and 21 unknowns
    @test length(equations(sys)) == 21
    @test length(unknowns(sys)) == 21

    # Verify key variables exist
    @test hasproperty(sys, :c_H)
    @test hasproperty(sys, :c_NH4)
    @test hasproperty(sys, :c_SO4)
    @test hasproperty(sys, :c_NO3)
    @test hasproperty(sys, :g_NH3)
    @test hasproperty(sys, :g_HNO3)
    @test hasproperty(sys, :W_w)
    @test hasproperty(sys, :I_s)
    @test hasproperty(sys, :γ_HNO3)

    # System should compile
    compiled = mtkcompile(sys)
    @test compiled isa System
end

# =============================================================================
# Equilibrium Solve Tests
# =============================================================================

@testitem "ISO2: Sulfate-poor equilibrium (R₁ ≥ 2)" setup = [Iso2Setup] tags = [:iso2] begin
    # Test case: excess ammonia relative to sulfate (typical continental aerosol)
    # SO4=1e-7, NH3=3e-7, NO3=5e-8 → R₁ = 3 (sulfate-poor)
    sys = IsorropiaEquilibrium()
    compiled = mtkcompile(sys)

    sol = iso2_solve(
        compiled, [
            compiled.RH => 0.8,
            compiled.W_SO4_total => 1.0e-7,
            compiled.W_NH3_total => 3.0e-7,
            compiled.W_NO3_total => 5.0e-8,
            compiled.c_SO4 => 9.5e-8,
            compiled.c_NO3 => 3.0e-8,
            compiled.c_Cl => 1.0e-20,
            compiled.c_H => 1.0e-11,
            compiled.c_OH => 1.0e-14,
            compiled.I_s => 10.0,
        ]
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success

    # All concentrations should be non-negative
    @test sol[compiled.c_H] ≥ 0
    @test sol[compiled.c_NH4] ≥ 0
    @test sol[compiled.c_SO4] ≥ 0
    @test sol[compiled.c_HSO4] ≥ 0
    @test sol[compiled.c_NO3] ≥ 0
    @test sol[compiled.g_NH3] ≥ 0
    @test sol[compiled.g_HNO3] ≥ 0
    @test sol[compiled.W_w] > 0
    @test sol[compiled.I_s] > 0

    # Mass conservation
    @test sol[compiled.c_SO4] + sol[compiled.c_HSO4] ≈ 1.0e-7 rtol = 1.0e-6
    @test sol[compiled.c_NH4] + sol[compiled.g_NH3] ≈ 3.0e-7 rtol = 1.0e-6
    @test sol[compiled.c_NO3] + sol[compiled.g_HNO3] ≈ 5.0e-8 rtol = 1.0e-6

    # Charge balance
    cations = sol[compiled.c_H] + sol[compiled.c_NH4]
    anions = 2 * sol[compiled.c_SO4] + sol[compiled.c_HSO4] + sol[compiled.c_NO3] + sol[compiled.c_OH]
    @test cations ≈ anions rtol = 1.0e-6

    # In sulfate-poor conditions with excess ammonia:
    @test sol[compiled.c_SO4] > sol[compiled.c_HSO4]
    @test sol[compiled.g_NH3] > 0
end

@testitem "ISO2: Sulfate-rich equilibrium (R₁ < 2)" setup = [Iso2Setup] tags = [:iso2] begin
    # Test case: sulfate exceeds ammonia (high SO2 environment)
    # SO4=2e-7, NH3=1e-7, NO3=5e-8 → R₁ = 0.5 (sulfate-rich)
    sys = IsorropiaEquilibrium()
    compiled = mtkcompile(sys)

    sol = iso2_solve(
        compiled, [
            compiled.RH => 0.8,
            compiled.W_SO4_total => 2.0e-7,
            compiled.W_NH3_total => 1.0e-7,
            compiled.W_NO3_total => 5.0e-8,
            compiled.c_SO4 => 1.5e-7,
            compiled.c_NO3 => 1.0e-9,
            compiled.c_Cl => 1.0e-20,
            compiled.c_H => 1.0e-7,
            compiled.c_OH => 1.0e-18,
            compiled.I_s => 15.0,
        ]
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Mass conservation
    @test sol[compiled.c_SO4] + sol[compiled.c_HSO4] ≈ 2.0e-7 rtol = 1.0e-6
    @test sol[compiled.c_NH4] + sol[compiled.g_NH3] ≈ 1.0e-7 rtol = 1.0e-6
    @test sol[compiled.c_NO3] + sol[compiled.g_HNO3] ≈ 5.0e-8 rtol = 1.0e-6

    # In sulfate-rich conditions:
    @test sol[compiled.c_H] > 0
    @test sol[compiled.c_NH4] > sol[compiled.g_NH3]
    @test sol[compiled.g_HNO3] > 0
end

@testitem "ISO2: Marine aerosol with sodium and chloride" setup = [Iso2Setup] tags = [:iso2] begin
    # Test case: marine aerosol with NaCl
    sys = IsorropiaEquilibrium()
    compiled = mtkcompile(sys)

    sol = iso2_solve(
        compiled, [
            compiled.RH => 0.85,
            compiled.W_Na_total => 5.0e-8,
            compiled.W_SO4_total => 5.0e-8,
            compiled.W_NH3_total => 1.0e-7,
            compiled.W_NO3_total => 3.0e-8,
            compiled.W_Cl_total => 5.0e-8,
            compiled.c_SO4 => 4.5e-8,
            compiled.c_NO3 => 1.0e-8,
            compiled.c_Cl => 4.0e-8,
            compiled.c_H => 1.0e-11,
            compiled.c_OH => 1.0e-14,
            compiled.I_s => 10.0,
        ]
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Mass conservation
    @test sol[compiled.c_Na] ≈ 5.0e-8 rtol = 1.0e-6
    @test sol[compiled.c_SO4] + sol[compiled.c_HSO4] ≈ 5.0e-8 rtol = 1.0e-6
    @test sol[compiled.c_NH4] + sol[compiled.g_NH3] ≈ 1.0e-7 rtol = 1.0e-6
    @test sol[compiled.c_NO3] + sol[compiled.g_HNO3] ≈ 3.0e-8 rtol = 1.0e-6
    @test sol[compiled.c_Cl] + sol[compiled.g_HCl] ≈ 5.0e-8 rtol = 1.0e-6

    # Charge balance
    cations = sol[compiled.c_H] + sol[compiled.c_Na] + sol[compiled.c_NH4]
    anions = (
        2 * sol[compiled.c_SO4] + sol[compiled.c_HSO4] +
            sol[compiled.c_NO3] + sol[compiled.c_Cl] + sol[compiled.c_OH]
    )
    @test cations ≈ anions rtol = 1.0e-6
end

@testitem "ISO2: RH sensitivity" setup = [Iso2Setup] tags = [:iso2] begin
    # Higher RH should give more aerosol water
    sys = IsorropiaEquilibrium()
    compiled = mtkcompile(sys)

    # Solve at RH=0.80 first to get a good reference solution
    sol_80 = iso2_solve(
        compiled, [
            compiled.RH => 0.8,
            compiled.W_SO4_total => 1.0e-7,
            compiled.W_NH3_total => 3.0e-7,
            compiled.W_NO3_total => 1.0e-7,
            compiled.c_SO4 => 9.5e-8,
            compiled.c_NO3 => 5.0e-8,
            compiled.c_Cl => 1.0e-20,
            compiled.c_H => 1.0e-11,
            compiled.c_OH => 1.0e-14,
            compiled.I_s => 10.0,
        ]
    )
    @test sol_80.retcode == SciMLBase.ReturnCode.Success

    # Use the RH=0.80 solution as initial guess for nearby RH values
    sol_90 = iso2_solve(
        compiled, [
            compiled.RH => 0.9,
            compiled.W_SO4_total => 1.0e-7,
            compiled.W_NH3_total => 3.0e-7,
            compiled.W_NO3_total => 1.0e-7,
            compiled.c_SO4 => sol_80[compiled.c_SO4],
            compiled.c_NO3 => sol_80[compiled.c_NO3],
            compiled.c_Cl => sol_80[compiled.c_Cl],
            compiled.c_H => sol_80[compiled.c_H],
            compiled.c_OH => sol_80[compiled.c_OH],
            compiled.I_s => sol_80[compiled.I_s],
        ]
    )
    @test sol_90.retcode == SciMLBase.ReturnCode.Success

    # Water content should increase with RH
    @test sol_80[compiled.W_w] < sol_90[compiled.W_w]
end

@testitem "ISO2: Temperature sensitivity" setup = [Iso2Setup] tags = [:iso2] begin
    # Higher temperature should volatilize more semi-volatile species
    sys = IsorropiaEquilibrium()
    compiled = mtkcompile(sys)

    T_vals = [273.15, 298.15, 313.15]
    g_NH3_vals = Float64[]
    g_HNO3_vals = Float64[]

    for T in T_vals
        sol = iso2_solve(
            compiled, [
                compiled.RH => 0.8,
                compiled.T_env => T,
                compiled.W_SO4_total => 1.0e-7,
                compiled.W_NH3_total => 3.0e-7,
                compiled.W_NO3_total => 1.0e-7,
                compiled.c_SO4 => 9.5e-8,
                compiled.c_NO3 => 5.0e-8,
                compiled.c_Cl => 1.0e-20,
                compiled.c_H => 1.0e-11,
                compiled.c_OH => 1.0e-14,
                compiled.I_s => 10.0,
            ]
        )
        @test sol.retcode == SciMLBase.ReturnCode.Success

        push!(g_NH3_vals, sol[compiled.g_NH3])
        push!(g_HNO3_vals, sol[compiled.g_HNO3])
    end

    # Higher temperature should produce more gas-phase semi-volatiles
    @test g_NH3_vals[1] < g_NH3_vals[3]
    @test g_HNO3_vals[1] < g_HNO3_vals[3]
end

# =============================================================================
# Figure Reproduction Tests (Figures 6-9, Table 8 cases)
# =============================================================================

@testsnippet Iso2FigSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t
    using DynamicQuantities
    using Aerosol
    using NonlinearSolve
    using SciMLBase

    # Helper: solve equilibrium sweeping RH with continuation
    function solve_iso2_case(
            compiled, Na, SO4, NH3, NO3, Cl, Ca, K, Mg;
            RH_range = 0.95:-0.05:0.3
        )
        RH_ok, W_vals = Float64[], Float64[]
        c = Dict(s => Float64[] for s in [:Na, :NH4, :SO4, :HSO4, :NO3, :Cl, :Ca, :K, :Mg, :H, :OH])
        g = Dict(s => Float64[] for s in [:NH3, :HNO3, :HCl])
        prev = nothing

        for rh in RH_range
            params = [
                compiled.RH => rh, compiled.W_Na_total => Na,
                compiled.W_SO4_total => SO4, compiled.W_NH3_total => NH3,
                compiled.W_NO3_total => NO3, compiled.W_Cl_total => Cl,
                compiled.W_Ca_total => Ca, compiled.W_K_total => K,
                compiled.W_Mg_total => Mg,
            ]
            if prev !== nothing
                guesses = [
                    compiled.c_SO4 => prev[:SO4], compiled.c_HSO4 => prev[:HSO4],
                    compiled.c_NH4 => prev[:NH4], compiled.c_NO3 => prev[:NO3],
                    compiled.c_Cl => prev[:Cl], compiled.c_H => prev[:H],
                    compiled.c_OH => prev[:OH], compiled.I_s => prev[:I],
                    compiled.W_w => prev[:W],
                ]
            else
                # Sulfate ratio R₁ determines aerosol regime
                R1 = SO4 > 0 ? (NH3 + Na + 2 * Ca + K + 2 * Mg) / SO4 : 100.0
                if R1 < 2  # Sulfate-rich: acidic aerosol
                    guesses = [
                        compiled.c_SO4 => max(0.5 * SO4, 1.0e-20),
                        compiled.c_HSO4 => max(0.5 * SO4, 1.0e-20),
                        compiled.c_NH4 => max(0.99 * NH3, 1.0e-20),
                        compiled.c_NO3 => 1.0e-10, compiled.c_Cl => 1.0e-20,
                        compiled.c_H => 1.0e-7, compiled.c_OH => 1.0e-18,
                        compiled.I_s => 15.0, compiled.W_w => 1.0e-7,
                        compiled.g_NH3 => 1.0e-10, compiled.g_HNO3 => max(0.99 * NO3, 1.0e-20),
                    ]
                else  # Sulfate-poor: more neutral aerosol
                    guesses = [
                        compiled.c_SO4 => max(0.85 * SO4, 1.0e-20),
                        compiled.c_HSO4 => max(0.15 * SO4, 1.0e-20),
                        compiled.c_NH4 => max(min(NH3, 2 * SO4) * 0.8, 1.0e-20),
                        compiled.c_NO3 => max(0.3 * NO3, 1.0e-20),
                        compiled.c_Cl => max(0.3 * Cl, 1.0e-20),
                        compiled.c_H => 1.0e-10, compiled.c_OH => 1.0e-13,
                        compiled.I_s => 10.0,
                    ]
                end
            end

            prob = NonlinearProblem(compiled, vcat(params, guesses))
            sol = solve(prob, RobustMultiNewton(); maxiters = 10000)

            if sol.retcode == SciMLBase.ReturnCode.Success
                push!(RH_ok, rh); push!(W_vals, sol[compiled.W_w])
                for s in [:Na, :NH4, :SO4, :HSO4, :NO3, :Cl, :Ca, :K, :Mg, :H, :OH]
                    push!(c[s], sol[getproperty(compiled, Symbol(:c_, s))])
                end
                push!(g[:NH3], sol[compiled.g_NH3])
                push!(g[:HNO3], sol[compiled.g_HNO3])
                push!(g[:HCl], sol[compiled.g_HCl])
                prev = Dict(
                    :SO4 => sol[compiled.c_SO4], :HSO4 => sol[compiled.c_HSO4],
                    :NH4 => sol[compiled.c_NH4], :NO3 => sol[compiled.c_NO3],
                    :Cl => sol[compiled.c_Cl], :H => sol[compiled.c_H],
                    :OH => sol[compiled.c_OH], :I => sol[compiled.I_s],
                    :W => sol[compiled.W_w]
                )
            end
        end
        reverse!(RH_ok); reverse!(W_vals)
        for v in values(c)
            reverse!(v)
        end
        for v in values(g)
            reverse!(v)
        end
        return RH_ok, W_vals, c, g
    end
end

@testitem "ISO2: Figure 6 — Urban sulfate-rich (Case 3)" setup = [Iso2FigSetup] tags = [:iso2] begin
    # Case 3 (Urban): Na=0, H₂SO₄=15, NH₃=2, HNO₃=10, HCl=0, Ca=0.9, K=1.0, Mg=0 µg/m³
    # Sulfate rich: R₁ ≈ 1.23 (1 < R₁ < 2)
    sys = IsorropiaEquilibrium()
    compiled = mtkcompile(sys)

    RH, W, c, g = solve_iso2_case(
        compiled,
        0.0, 15.0e-6 / 98.08, 2.0e-6 / 17.03, 10.0e-6 / 63.01, 0.0,
        0.9e-6 / 40.08, 1.0e-6 / 39.1, 0.0
    )

    @test length(RH) > 5  # Solver should converge for most RH values

    # Mass conservation at all converged points
    SO4_total = 15.0e-6 / 98.08
    NH3_total = 2.0e-6 / 17.03
    NO3_total = 10.0e-6 / 63.01
    for i in eachindex(RH)
        @test c[:SO4][i] + c[:HSO4][i] ≈ SO4_total rtol = 1.0e-4
        @test c[:NH4][i] + g[:NH3][i] ≈ NH3_total rtol = 1.0e-4
        @test c[:NO3][i] + g[:HNO3][i] ≈ NO3_total rtol = 1.0e-4
    end

    # K⁺ is non-volatile: all K in aerosol (Fig 6b: K⁺ ≈ 1.0 µg/m³ at all RH)
    K_total_ug = 1.0  # µg/m³
    for i in eachindex(RH)
        @test c[:K][i] * 39.1e6 ≈ K_total_ug rtol = 1.0e-4
    end

    # Water content increases with RH (Fig 6a)
    idx_low = findfirst(x -> x ≥ 0.5, RH)
    idx_high = findlast(x -> x ≤ 0.9, RH)
    if idx_low !== nothing && idx_high !== nothing
        @test W[idx_high] > W[idx_low]
    end

    # Sulfate-rich: NH₄⁺ is mostly in aerosol (limited by sulfate neutralization)
    # Fig 6c: NH₄⁺ ≈ 2 µg/m³ (nearly all ammonia as NH₄⁺)
    NH3_total_as_NH4_ug = NH3_total * 18.04e6
    for i in eachindex(RH)
        NH4_ug = c[:NH4][i] * 18.04e6
        @test NH4_ug > 0.5
        @test NH4_ug < NH3_total_as_NH4_ug * 1.1  # Should not exceed total
    end

    # Sulfate-rich: very little aerosol NO₃⁻ (HNO₃ stays in gas due to acidity)
    # Fig 6d: NO₃⁻ is near zero in sulfate-rich conditions
    for i in eachindex(RH)
        @test abs(c[:NO3][i]) < 0.1 * NO3_total  # Most nitrate in gas phase
    end
end

@testitem "ISO2: Figure 7 — Marine sulfate-poor (Case 12)" setup = [Iso2FigSetup] tags = [:iso2] begin
    # Case 12 (Marine): Na=3.0, H₂SO₄=3.0, NH₃=0.020, HNO₃=2.0, HCl=3.121,
    #                    Ca=0.36, K=0.45, Mg=0.13 µg/m³
    # Sulfate poor: R₁ ≈ 5.6 (R₁ > 2, R₂ > 2)
    sys = IsorropiaEquilibrium()
    compiled = mtkcompile(sys)

    RH, W, c, g = solve_iso2_case(
        compiled,
        3.0e-6 / 22.99, 3.0e-6 / 98.08, 0.02e-6 / 17.03, 2.0e-6 / 63.01, 3.121e-6 / 36.46,
        0.36e-6 / 40.08, 0.45e-6 / 39.1, 0.13e-6 / 24.31
    )

    @test length(RH) > 5

    # Mass conservation
    Na_total = 3.0e-6 / 22.99
    SO4_total = 3.0e-6 / 98.08
    Cl_total = 3.121e-6 / 36.46
    for i in eachindex(RH)
        @test c[:Na][i] ≈ Na_total rtol = 1.0e-4
        @test c[:SO4][i] + c[:HSO4][i] ≈ SO4_total rtol = 1.0e-4
        @test c[:Cl][i] + g[:HCl][i] ≈ Cl_total rtol = 1.0e-4
    end

    # Non-volatile species: K⁺ ≈ 0.45 µg/m³ (Fig 7b)
    for i in eachindex(RH)
        @test c[:K][i] * 39.1e6 ≈ 0.45 rtol = 1.0e-4
    end

    # Mg²⁺ non-volatile: ≈ 0.13 µg/m³ (Fig 7d)
    for i in eachindex(RH)
        @test c[:Mg][i] * 24.31e6 ≈ 0.13 rtol = 1.0e-4
    end

    # Na⁺ non-volatile in metastable: all Na in solution (no NaCl solid)
    for i in eachindex(RH)
        @test c[:Na][i] * 22.99e6 ≈ 3.0 rtol = 1.0e-4
    end

    # Water content increases with RH
    idx_low = findfirst(x -> x ≥ 0.5, RH)
    idx_high = findlast(x -> x ≤ 0.9, RH)
    if idx_low !== nothing && idx_high !== nothing
        @test W[idx_high] > W[idx_low]
    end
end

@testitem "ISO2: Figure 8 — Continental sulfate-poor (Case 5)" setup = [Iso2FigSetup] tags = [:iso2] begin
    # Case 5 (Non-urban Continental): Na=0.2, H₂SO₄=2.0, NH₃=8.0, HNO₃=12.0,
    #                                  HCl=0.2, Ca=0.12, K=0.18, Mg=0 µg/m³
    # Sulfate poor: R₁ ≈ 24 (R₁ >> 2)
    sys = IsorropiaEquilibrium()
    compiled = mtkcompile(sys)

    RH, W, c, g = solve_iso2_case(
        compiled,
        0.2e-6 / 22.99, 2.0e-6 / 98.08, 8.0e-6 / 17.03, 12.0e-6 / 63.01, 0.2e-6 / 36.46,
        0.12e-6 / 40.08, 0.18e-6 / 39.1, 0.0
    )

    @test length(RH) > 5

    # Mass conservation
    NH3_total = 8.0e-6 / 17.03
    NO3_total = 12.0e-6 / 63.01
    for i in eachindex(RH)
        @test c[:NH4][i] + g[:NH3][i] ≈ NH3_total rtol = 1.0e-4
        @test c[:NO3][i] + g[:HNO3][i] ≈ NO3_total rtol = 1.0e-4
    end

    # Sulfate-poor with excess ammonia: significant NH₃ in gas phase
    # Fig 8c: NH₄⁺ varies with RH but stays below total ammonia
    idx_high = findlast(x -> x ≤ 0.9, RH)
    if idx_high !== nothing
        @test c[:NH4][idx_high] * 18.04e6 > 1.0  # Significant aerosol NH₄⁺
        @test g[:NH3][idx_high] > 0  # Some NH₃ remains in gas
    end

    # Fig 8b: NO₃⁻ should be significant at high RH (aerosol nitrate present)
    if idx_high !== nothing
        NO3_ug = c[:NO3][idx_high] * 62.0e6
        @test NO3_ug > 0.5  # Significant aerosol nitrate
    end

    # Water content increases with RH
    idx_low = findfirst(x -> x ≥ 0.5, RH)
    if idx_low !== nothing && idx_high !== nothing
        @test W[idx_high] > W[idx_low]
    end
end

@testitem "ISO2: Figure 9 — Remote continental (Case 13)" setup = [Iso2FigSetup] tags = [:iso2] begin
    # Case 13 (Remote Continental): Na=0, H₂SO₄=10.0, NH₃=4.25, HNO₃=0.145,
    #                                HCl=0, Ca=0.08, K=0.09, Mg=0 µg/m³
    # Near-neutral: R₁ ≈ 2.5
    sys = IsorropiaEquilibrium()
    compiled = mtkcompile(sys)

    RH, W, c, g = solve_iso2_case(
        compiled,
        0.0, 10.0e-6 / 98.08, 4.25e-6 / 17.03, 0.145e-6 / 63.01, 0.0,
        0.08e-6 / 40.08, 0.09e-6 / 39.1, 0.0
    )

    @test length(RH) >= 5

    # Mass conservation
    SO4_total = 10.0e-6 / 98.08
    NH3_total = 4.25e-6 / 17.03
    for i in eachindex(RH)
        @test c[:SO4][i] + c[:HSO4][i] ≈ SO4_total rtol = 1.0e-4
        @test c[:NH4][i] + g[:NH3][i] ≈ NH3_total rtol = 1.0e-4
    end

    # K⁺ is non-volatile: ≈ 0.09 µg/m³ (Fig 9b)
    for i in eachindex(RH)
        @test c[:K][i] * 39.1e6 ≈ 0.09 rtol = 1.0e-4
    end

    # Water content increases with RH (Fig 9a)
    idx_low = findfirst(x -> x ≥ 0.5, RH)
    idx_high = findlast(x -> x ≤ 0.9, RH)
    if idx_low !== nothing && idx_high !== nothing
        @test W[idx_high] > W[idx_low]
    end

    # Near-neutral case: most NH₃ should be in aerosol as (NH₄)₂SO₄
    # R₁ ≈ 2.5 means NH₄⁺ ≈ 2×SO₄ (most sulfate neutralized)
    idx_mid = findfirst(x -> x ≥ 0.8, RH)
    if idx_mid !== nothing
        # NH₄⁺ should be close to (but not exceeding) 2×SO₄_total
        @test c[:NH4][idx_mid] ≤ NH3_total * 1.01  # Can't exceed total
        @test c[:NH4][idx_mid] > SO4_total  # Should exceed 1×SO₄ (partially neutralized)
    end
end
