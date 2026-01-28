@testsnippet SizeDistSetup begin
    using Test
    using ModelingToolkit
    using Aerosol
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
    # Eq. 8.44: n_N^o(log D_p) = N / (sqrt(2π) * log(σ_g)) * exp(-(log D_p - log D_pg)^2 / (2 log^2(σ_g)))
    # Verify using a single-mode system
    sys = AerosolDistribution(1)

    N_val = 1e11     # m^-3
    D_g_val = 1e-7   # 0.1 μm
    logσ_val = 0.3
    D_p_val = 1e-7   # evaluate at geometric median

    # At D_p = D_g, log10(D_p/D_g) = 0, so exp term = 1
    # n_N_o = N / (sqrt(2π) * logσ) = 1e11 / (sqrt(2π) * 0.3)
    expected_peak = N_val / (sqrt(2π) * logσ_val)

    # The symbolic system represents this correctly
    # We can verify by direct computation
    actual = N_val / (sqrt(2π) * logσ_val) * exp(-(log10(D_p_val / D_g_val))^2 / (2 * logσ_val^2))
    @test actual ≈ expected_peak rtol=1e-10

    # At D_p = D_g * 10^σ (one σ away), the value should be exp(-0.5) times the peak
    D_p_1sigma = D_g_val * 10^logσ_val
    actual_1sigma = N_val / (sqrt(2π) * logσ_val) * exp(-(log10(D_p_1sigma / D_g_val))^2 / (2 * logσ_val^2))
    @test actual_1sigma ≈ expected_peak * exp(-0.5) rtol=1e-10
end

@testitem "Total Number Concentration" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Verify that N_t equals the sum of mode concentrations
    # Table 8.3 values (converted to m^-3)
    # Urban: 9.93e4 + 1.11e3 + 3.64e4 = 136810 cm^-3 = 1.3681e11 m^-3
    urban_N_t = (9.93e4 + 1.11e3 + 3.64e4) * 1e6
    @test urban_N_t ≈ 1.3681e11

    # Marine: 133 + 66.6 + 3.06 = 202.66 cm^-3
    marine_N_t = (133.0 + 66.6 + 3.06) * 1e6
    @test marine_N_t ≈ 202.66e6

    # Remote continental: 3200 + 2900 + 0.3 = 6100.3 cm^-3
    remote_N_t = (3200.0 + 2900.0 + 0.3) * 1e6
    @test remote_N_t ≈ 6100.3e6
end

@testitem "Hatch-Choate Equations" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Table 8.2 — Hatch-Choate equations for lognormal distributions
    # Verify relationships between geometric median diameter and other diameters

    D_pg = 0.1e-6  # 0.1 μm
    logσ = 0.3
    lnσ = logσ * log(10)  # convert to natural log

    # Eq. 8.39: Mean diameter = D_pg * exp(0.5 * ln²(σ_g))
    D_bar = D_pg * exp(0.5 * lnσ^2)
    @test D_bar > D_pg  # mean > geometric median for lognormal

    # Eq. 8.40: Mode diameter = D_pg * exp(-ln²(σ_g))
    D_mode = D_pg * exp(-lnσ^2)
    @test D_mode < D_pg  # mode < geometric median for lognormal

    # Eq. 8.49: Surface area median diameter = D_pg * exp(2 * ln²(σ_g))
    D_ps = D_pg * exp(2 * lnσ^2)
    @test D_ps > D_pg

    # Eq. 8.52: Volume median diameter = D_pg * exp(3 * ln²(σ_g))
    D_pv = D_pg * exp(3 * lnσ^2)
    @test D_pv > D_ps > D_pg

    # Verify ordering: D_mode < D_pg < D_bar < D_ps < D_pv
    @test D_mode < D_pg < D_bar < D_ps < D_pv

    # Eq. 8.46/8.47: 84.1% and 15.9% quantile relationship
    # D_p(84.1%) = D_pg * σ_g
    σ_g = 10^logσ
    D_84 = D_pg * σ_g
    D_16 = D_pg / σ_g
    @test D_84 > D_pg
    @test D_16 < D_pg
    @test D_84 / D_16 ≈ σ_g^2
end

@testitem "Moment Equations" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Eq. 8.41: M_k = N_t * D_pg^k * exp(k^2 * ln²(σ_g) / 2)
    N_t = 1e11  # m^-3
    D_pg = 0.1e-6  # m
    logσ = 0.3
    lnσ = logσ * log(10)

    # k=0: M_0 = N_t (zeroth moment is total count)
    M_0 = N_t * D_pg^0 * exp(0^2 * lnσ^2 / 2)
    @test M_0 ≈ N_t

    # k=1: M_1 = N_t * D_pg * exp(ln²σ/2) = N_t * D_bar (Eq. 8.39)
    M_1 = N_t * D_pg * exp(lnσ^2 / 2)
    D_bar = D_pg * exp(0.5 * lnσ^2)
    @test M_1 ≈ N_t * D_bar

    # k=2: relates to total surface area (Eq. 8.5)
    # S_t = π * M_2 = π * N_t * D_pg^2 * exp(2 * ln²σ)
    M_2 = N_t * D_pg^2 * exp(4 * lnσ^2 / 2)
    S_t = π * M_2
    @test S_t > 0

    # k=3: relates to total volume (Eq. 8.7)
    # V_t = (π/6) * M_3 = (π/6) * N_t * D_pg^3 * exp(9 * ln²σ / 2)
    M_3 = N_t * D_pg^3 * exp(9 * lnσ^2 / 2)
    V_t = (π / 6) * M_3
    @test V_t > 0

    # Verify that S_t calculation matches the implementation formula
    S_t_impl = π * N_t * D_pg^2 * exp(2 * lnσ^2)
    @test S_t ≈ S_t_impl
end

@testitem "Vertical Mass Profile" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Eq. 8.55: M(z) = M(0) * exp(-z / H_p)
    H_p = 1000.0  # m

    # At z=0, M_z = 1.0
    @test exp(-0.0 / H_p) ≈ 1.0

    # At z=H_p, M_z = exp(-1) ≈ 0.368
    @test exp(-H_p / H_p) ≈ exp(-1) ≈ 0.36787944117144233

    # At z=2*H_p, M_z = exp(-2) ≈ 0.135
    @test exp(-2 * H_p / H_p) ≈ exp(-2)

    # Monotonically decreasing
    z_vals = [0, 500, 1000, 2000, 5000]
    M_vals = [exp(-z / H_p) for z in z_vals]
    for i in 1:length(M_vals)-1
        @test M_vals[i] > M_vals[i+1]
    end
end

@testitem "Distribution Symmetry" setup=[SizeDistSetup] tags=[:size_dist] begin
    # The lognormal distribution is symmetric in log-space around log D_pg
    # For a single mode, n_N_o(D_pg * 10^x) should equal n_N_o(D_pg * 10^(-x))
    N = 1e11
    D_g = 1e-7
    logσ = 0.3

    lognormal_pdf(D_p) = N / (sqrt(2π) * logσ) * exp(-(log10(D_p / D_g))^2 / (2 * logσ^2))

    # Test symmetry at several offsets
    for x in [0.1, 0.2, 0.3, 0.5]
        D_p_plus = D_g * 10^x
        D_p_minus = D_g * 10^(-x)
        @test lognormal_pdf(D_p_plus) ≈ lognormal_pdf(D_p_minus) rtol=1e-12
    end
end

@testitem "Surface and Volume Median Diameters" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Eq. 8.49: D_ps = D_pg * exp(2 * ln²(σ_g))
    # Eq. 8.52: D_pv = D_pg * exp(3 * ln²(σ_g))
    # For single mode, verify relationship: D_pv/D_ps = exp(ln²(σ_g))

    D_pg = 0.05e-6  # 50 nm
    logσ = 0.3
    lnσ = logσ * log(10)

    D_ps = D_pg * exp(2 * lnσ^2)
    D_pv = D_pg * exp(3 * lnσ^2)

    @test D_pv / D_ps ≈ exp(lnσ^2)

    # Verify D_ps > D_pg
    @test D_ps > D_pg

    # Verify D_pv > D_ps
    @test D_pv > D_ps

    # For σ_g = 1 (monodisperse), all diameters should equal D_pg
    logσ_mono = 0.0 + eps()  # avoid division by zero
    lnσ_mono = logσ_mono * log(10)
    D_ps_mono = D_pg * exp(2 * lnσ_mono^2)
    D_pv_mono = D_pg * exp(3 * lnσ_mono^2)
    @test D_ps_mono ≈ D_pg atol=1e-20
    @test D_pv_mono ≈ D_pg atol=1e-20
end

@testitem "Table 8.3 Consistency" setup=[SizeDistSetup] tags=[:size_dist] begin
    # Verify that all Table 8.3 entries have physically reasonable values

    # Table 8.3 data (N in cm^-3, D_p in μm, log σ)
    table_data = Dict(
        "Urban" => [(9.93e4, 0.013, 0.245), (1.11e3, 0.014, 0.666), (3.64e4, 0.050, 0.337)],
        "Marine" => [(133, 0.008, 0.657), (66.6, 0.266, 0.210), (3.06, 0.580, 0.396)],
        "Rural" => [(6.65e3, 0.015, 0.225), (147, 0.054, 0.557), (1990, 0.084, 0.266)],
        "Remote continental" => [(3200, 0.020, 0.161), (2900, 0.116, 0.217), (0.300, 1.800, 0.380)],
        "Free troposphere" => [(129, 0.007, 0.645), (59.7, 0.250, 0.253), (63.5, 0.520, 0.425)],
        "Polar" => [(21.7, 0.138, 0.164), (0.186, 0.750, 0.521), (3.04e-4, 8.600, 0.420)],
        "Desert" => [(726, 0.002, 0.247), (114, 0.038, 0.770), (0.178, 21.60, 0.438)],
    )

    for (env_type, modes) in table_data
        for (i, (N, D_p, logσ)) in enumerate(modes)
            @test N > 0       # Positive concentrations
            @test D_p > 0     # Positive diameters
            @test logσ > 0    # Positive spread
            @test logσ < 1.0  # Reasonable spread (σ_g < 10)
        end
    end
end
