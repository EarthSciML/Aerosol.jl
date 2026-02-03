@testsnippet SizeDistSetup begin
    using Test
    using ModelingToolkit
    using Symbolics
    using Aerosol

    """
    Helper to get the base name of a symbolic variable (part after last ₊),
    and strip the (t) suffix if present.
    """
    function get_base_name(sym)
        name = string(Symbolics.tosymbol(sym, escape = false))
        name = contains(name, "₊") ? split(name, "₊")[end] : name
        name = replace(name, r"\(t\)$" => "")
        return name
    end

    """
    Helper to evaluate a system equation's RHS by substituting parameter values.
    Returns the numerical result of the RHS for the equation whose LHS matches `var`.

    This function handles namespaced variables by matching the actual symbols
    in the equation RHS with the user-provided parameter values.
    """
    function eval_eq(sys, var, param_vals)
        var_base = get_base_name(var)

        # Find the equation for this variable
        matching_eqs = filter(equations(sys)) do eq
            lhs_base = get_base_name(eq.lhs)
            return lhs_base == var_base
        end
        eq = only(matching_eqs)

        # Get all variables that appear in the RHS
        rhs_vars = Symbolics.get_variables(eq.rhs)

        # Build substitution dict by matching RHS variables with:
        # 1. System defaults (for constants like π_c, ln10)
        # 2. User-provided values (matching by string representation)
        subs = Dict{Any,Any}()
        defs = ModelingToolkit.defaults(sys)

        for rv in rhs_vars
            rv_str = string(rv)

            # Check system defaults
            for (k, v) in defs
                k_str = string(k)
                # Match if rv equals k, or k is namespaced version of rv
                if rv_str == k_str || endswith(k_str, "₊" * rv_str)
                    subs[rv] = v
                end
            end

            # Check user-provided values (these override defaults)
            for (k, v) in param_vals
                k_str = string(k)
                # Match if rv equals k, or k is namespaced version of rv
                if rv_str == k_str || endswith(k_str, "₊" * rv_str)
                    subs[rv] = v
                end
            end
        end

        result = Symbolics.substitute(eq.rhs, subs)
        return Float64(result)
    end
end

@testitem "Structural Verification" setup=[SizeDistSetup] tags=[:size_dist] begin
    sys = AerosolDistribution(3)
    @test length(equations(sys)) == 10
    @test length(unknowns(sys)) == 10

    # Verify key unknowns exist
    unknown_names = [string(v) for v in unknowns(sys)]
    @test "n_N_o(t)" in unknown_names
    @test "n_S_o(t)" in unknown_names
    @test "n_V_o(t)" in unknown_names
    @test "N_t(t)" in unknown_names
    @test "S_t(t)" in unknown_names
    @test "V_t(t)" in unknown_names
    @test "D_s(t)" in unknown_names
    @test "D_v(t)" in unknown_names
    @test "D_bar(t)" in unknown_names
    @test "M_z(t)" in unknown_names

    # Verify key parameters exist
    param_names = [string(p) for p in parameters(sys)]
    @test any(contains("N"), param_names)
    @test any(contains("D_g"), param_names)
    @test any(contains("logσ"), param_names)
    @test any(contains("D_p"), param_names)
    @test any(contains("z"), param_names)
    @test any(contains("H_p"), param_names)

    # Test different number of modes
    sys2 = AerosolDistribution(2)
    @test length(equations(sys2)) == 10

    sys5 = AerosolDistribution(5)
    @test length(equations(sys5)) == 10
end

@testitem "Predefined Distributions" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Test that all 7 predefined distributions can be created
    distributions = [
        UrbanAerosol,
        MarineAerosol,
        RuralAerosol,
        RemoteContinentalAerosol,
        FreeTroposphereAerosol,
        PolarAerosol,
        DesertAerosol,
    ]

    for dist_fn in distributions
        sys = dist_fn()
        @test length(equations(sys)) == 10
        defs = ModelingToolkit.defaults(sys)
        @test length(defs) > 0
    end

    # Verify specific default values (Table 8.3 — Urban)
    urban = UrbanAerosol()
    defs = ModelingToolkit.defaults(urban)

    # Find the N parameter defaults (converted from cm^-3 to m^-3)
    n_keys = filter(k -> contains(string(k), "N["), collect(keys(defs)))
    n_vals = sort([defs[k] for k in n_keys], rev=true)
    @test n_vals[1] ≈ 9.93e4 * 1e6  # Mode 1 N
    @test n_vals[2] ≈ 3.64e4 * 1e6  # Mode 3 N
    @test n_vals[3] ≈ 1.11e3 * 1e6  # Mode 2 N
end

@testitem "Lognormal Number Distribution" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Eq. 8.54: Verify the system's n_N_o equation evaluates correctly
    sys = AerosolDistribution(1)

    N_val = 1e11     # m^-3
    D_g_val = 1e-7   # 0.1 μm
    logσ_val = 0.3
    base_params = Dict(sys.N[1] => N_val, sys.D_g[1] => D_g_val, sys.logσ[1] => logσ_val)

    # At D_p = D_g, log10(D_p/D_g) = 0, so exp term = 1
    # n_N_o = N / (sqrt(2π) * logσ)
    expected_peak = N_val / (sqrt(2π) * logσ_val)
    result = eval_eq(sys, sys.n_N_o, merge(base_params, Dict(sys.D_p => D_g_val)))
    @test result ≈ expected_peak rtol=1e-10

    # At D_p = D_g * 10^σ (one σ away), the value should be exp(-0.5) times the peak
    D_p_1sigma = D_g_val * 10^logσ_val
    result_1sig = eval_eq(sys, sys.n_N_o, merge(base_params, Dict(sys.D_p => D_p_1sigma)))
    @test result_1sig ≈ expected_peak * exp(-0.5) rtol=1e-10
end

@testitem "Total Number Concentration" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Verify the system's N_t equation sums mode concentrations correctly
    sys = AerosolDistribution(3)

    # Use Urban Table 8.3 values (converted to m^-3)
    N_vals = [9.93e4 * 1e6, 1.11e3 * 1e6, 3.64e4 * 1e6]
    params = Dict(
        sys.N[1] => N_vals[1], sys.N[2] => N_vals[2], sys.N[3] => N_vals[3],
        sys.D_g[1] => 0.013e-6, sys.D_g[2] => 0.014e-6, sys.D_g[3] => 0.050e-6,
        sys.logσ[1] => 0.245, sys.logσ[2] => 0.666, sys.logσ[3] => 0.337,
        sys.D_p => 1e-7,
    )

    result = eval_eq(sys, sys.N_t, params)
    @test result ≈ sum(N_vals) rtol=1e-10

    # Marine: 133 + 66.6 + 3.1 = 202.7 cm^-3 (Table 8.3)
    marine = MarineAerosol()
    marine_params = Dict(marine.D_p => 1e-7)
    marine_Nt = eval_eq(marine, marine.N_t, marine_params)
    @test marine_Nt ≈ (133.0 + 66.6 + 3.1) * 1e6 rtol=1e-10

    # Remote continental: 3200 + 2900 + 0.3 = 6100.3 cm^-3
    remote = RemoteContinentalAerosol()
    remote_params = Dict(remote.D_p => 1e-7)
    remote_Nt = eval_eq(remote, remote.N_t, remote_params)
    @test remote_Nt ≈ (3200.0 + 2900.0 + 0.3) * 1e6 rtol=1e-10
end

@testitem "Hatch-Choate Equations" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Verify the system's D_s, D_v, D_bar equations implement Table 8.2 correctly
    sys = AerosolDistribution(1)

    D_pg = 0.1e-6  # 0.1 μm
    logσ = 0.3
    lnσ = logσ * log(10)  # convert to natural log
    params = Dict(
        sys.N[1] => 1e11, sys.D_g[1] => D_pg, sys.logσ[1] => logσ,
        sys.D_p => D_pg,
    )

    # Eq. 8.39: Mean diameter = D_pg * exp(0.5 * ln²(σ_g))
    D_bar_result = eval_eq(sys, sys.D_bar, params)
    @test D_bar_result ≈ D_pg * exp(0.5 * lnσ^2) rtol=1e-10
    @test D_bar_result > D_pg  # mean > geometric median for lognormal

    # Eq. 8.49: Surface area median diameter = D_pg * exp(2 * ln²(σ_g))
    D_s_result = eval_eq(sys, sys.D_s, params)
    @test D_s_result ≈ D_pg * exp(2 * lnσ^2) rtol=1e-10
    @test D_s_result > D_pg

    # Eq. 8.52: Volume median diameter = D_pg * exp(3 * ln²(σ_g))
    D_v_result = eval_eq(sys, sys.D_v, params)
    @test D_v_result ≈ D_pg * exp(3 * lnσ^2) rtol=1e-10
    @test D_v_result > D_s_result > D_pg

    # Verify ordering: D_pg < D_bar < D_s < D_v
    @test D_pg < D_bar_result < D_s_result < D_v_result

    # Verify D_v / D_s = exp(ln²σ)
    @test D_v_result / D_s_result ≈ exp(lnσ^2) rtol=1e-10
end

@testitem "Moment Equations" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Verify the system's S_t and V_t equations implement moment-based formulas
    # Eq. 8.41: M_k = N_t * D_pg^k * exp(k^2 * ln²(σ_g) / 2)
    sys = AerosolDistribution(1)

    N_val = 1e11  # m^-3
    D_pg = 0.1e-6  # m
    logσ = 0.3
    lnσ = logσ * log(10)
    params = Dict(
        sys.N[1] => N_val, sys.D_g[1] => D_pg, sys.logσ[1] => logσ,
        sys.D_p => D_pg,
    )

    # Total surface area — Eq. 8.5 with Eq. 8.41 (k=2 moment)
    # S_t = π * N * D_pg^2 * exp(2 * ln²σ)
    S_t_result = eval_eq(sys, sys.S_t, params)
    expected_S_t = π * N_val * D_pg^2 * exp(2 * lnσ^2)
    @test S_t_result ≈ expected_S_t rtol=1e-10
    @test S_t_result > 0

    # Total volume — Eq. 8.7 with Eq. 8.41 (k=3 moment)
    # V_t = (π/6) * N * D_pg^3 * exp(4.5 * ln²σ)
    V_t_result = eval_eq(sys, sys.V_t, params)
    expected_V_t = (π / 6) * N_val * D_pg^3 * exp(4.5 * lnσ^2)
    @test V_t_result ≈ expected_V_t rtol=1e-10
    @test V_t_result > 0

    # Verify N_t for single mode equals N
    N_t_result = eval_eq(sys, sys.N_t, params)
    @test N_t_result ≈ N_val rtol=1e-10
end

@testitem "Vertical Mass Profile" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Eq. 8.55: Verify the system's M_z equation
    sys = AerosolDistribution(1)
    base_params = Dict(
        sys.N[1] => 1e11, sys.D_g[1] => 1e-7, sys.logσ[1] => 0.3,
        sys.D_p => 1e-7,
    )

    # At z=0, M_z = 1.0
    result_z0 = eval_eq(sys, sys.M_z, merge(base_params, Dict(sys.z => 0.0, sys.H_p => 1000.0)))
    @test result_z0 ≈ 1.0 rtol=1e-10

    # At z=H_p, M_z = exp(-1) ≈ 0.368
    H_p = 1000.0
    result_zHp = eval_eq(sys, sys.M_z, merge(base_params, Dict(sys.z => H_p, sys.H_p => H_p)))
    @test result_zHp ≈ exp(-1) rtol=1e-10

    # At z=2*H_p, M_z = exp(-2) ≈ 0.135
    result_z2Hp = eval_eq(sys, sys.M_z, merge(base_params, Dict(sys.z => 2 * H_p, sys.H_p => H_p)))
    @test result_z2Hp ≈ exp(-2) rtol=1e-10

    # Monotonically decreasing
    z_vals = [0.0, 500.0, 1000.0, 2000.0, 5000.0]
    M_vals = [eval_eq(sys, sys.M_z, merge(base_params, Dict(sys.z => z, sys.H_p => H_p))) for z in z_vals]
    for i in 1:length(M_vals)-1
        @test M_vals[i] > M_vals[i+1]
    end
end

@testitem "Distribution Symmetry" setup=[SizeDistSetup] tags=[:size_dist] begin
    # The lognormal distribution is symmetric in log-space around log D_pg
    # For a single mode, n_N_o(D_pg * 10^x) should equal n_N_o(D_pg * 10^(-x))
    sys = AerosolDistribution(1)

    N_val = 1e11
    D_g_val = 1e-7
    logσ_val = 0.3
    base_params = Dict(sys.N[1] => N_val, sys.D_g[1] => D_g_val, sys.logσ[1] => logσ_val)

    # Test symmetry at several offsets using the system's n_N_o equation
    for x in [0.1, 0.2, 0.3, 0.5]
        D_p_plus = D_g_val * 10^x
        D_p_minus = D_g_val * 10^(-x)
        val_plus = eval_eq(sys, sys.n_N_o, merge(base_params, Dict(sys.D_p => D_p_plus)))
        val_minus = eval_eq(sys, sys.n_N_o, merge(base_params, Dict(sys.D_p => D_p_minus)))
        @test val_plus ≈ val_minus rtol=1e-12
    end
end

@testitem "Surface and Volume Median Diameters" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Verify the system's D_s and D_v equations for various parameter values
    sys = AerosolDistribution(1)

    D_pg = 0.05e-6  # 50 nm
    logσ = 0.3
    lnσ = logσ * log(10)
    params = Dict(
        sys.N[1] => 1e11, sys.D_g[1] => D_pg, sys.logσ[1] => logσ,
        sys.D_p => D_pg,
    )

    D_s_result = eval_eq(sys, sys.D_s, params)
    D_v_result = eval_eq(sys, sys.D_v, params)

    # Eq. 8.49 / 8.52: D_pv / D_ps = exp(ln²σ)
    @test D_v_result / D_s_result ≈ exp(lnσ^2) rtol=1e-10

    # Verify D_ps > D_pg
    @test D_s_result > D_pg

    # Verify D_pv > D_ps
    @test D_v_result > D_s_result

    # For near-monodisperse (logσ → 0), all diameters should approach D_pg
    logσ_mono = 1e-6
    params_mono = Dict(
        sys.N[1] => 1e11, sys.D_g[1] => D_pg, sys.logσ[1] => logσ_mono,
        sys.D_p => D_pg,
    )
    D_s_mono = eval_eq(sys, sys.D_s, params_mono)
    D_v_mono = eval_eq(sys, sys.D_v, params_mono)
    @test D_s_mono ≈ D_pg rtol=1e-6
    @test D_v_mono ≈ D_pg rtol=1e-6
end

@testitem "Table 8.3 Consistency" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Verify that all predefined distributions have physically reasonable defaults
    distributions = Dict(
        "Urban" => UrbanAerosol,
        "Marine" => MarineAerosol,
        "Rural" => RuralAerosol,
        "Remote continental" => RemoteContinentalAerosol,
        "Free troposphere" => FreeTroposphereAerosol,
        "Polar" => PolarAerosol,
        "Desert" => DesertAerosol,
    )

    for (env_type, dist_fn) in distributions
        sys = dist_fn()
        defs = ModelingToolkit.defaults(sys)

        # Check that N, D_g, logσ defaults exist and are physical
        for i in 1:3
            n_keys = filter(k -> contains(string(k), "N[$i]"), collect(keys(defs)))
            d_keys = filter(k -> contains(string(k), "D_g[$i]"), collect(keys(defs)))
            s_keys = filter(k -> contains(string(k), "logσ[$i]"), collect(keys(defs)))
            @test length(n_keys) == 1
            @test length(d_keys) == 1
            @test length(s_keys) == 1
            @test defs[n_keys[1]] > 0       # Positive concentrations
            @test defs[d_keys[1]] > 0       # Positive diameters
            @test defs[s_keys[1]] > 0       # Positive spread
            @test defs[s_keys[1]] < 1.0     # Reasonable spread (σ_g < 10)
        end

        # Verify N_t evaluates to the sum of mode N values
        params = Dict(sys.D_p => 1e-7)
        N_t_result = eval_eq(sys, sys.N_t, params)
        n_keys = filter(k -> contains(string(k), "N["), collect(keys(defs)))
        N_sum = sum(defs[k] for k in n_keys)
        @test N_t_result ≈ N_sum rtol=1e-10
    end
end
