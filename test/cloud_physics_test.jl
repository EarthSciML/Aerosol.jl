@testitem "WaterProperties structural" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = WaterProperties()
    @test length(equations(sys)) == 7
    @test length(unknowns(sys)) == 7

    # Check that key variable names exist
    var_names = string.(unknowns(sys))
    @test any(contains("T_C"), var_names)
    @test any(contains("c_pw"), var_names)
    @test any(contains("ΔH_v"), var_names)
    @test any(contains("ΔH_m"), var_names)
    @test any(contains("σ_w0"), var_names)
    @test any(contains("p_sat_water"), var_names)
    @test any(contains("p_sat_ice"), var_names)
end

@testitem "WaterProperties saturation vapor pressure" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = WaterProperties()
    ssys = mtkcompile(sys)

    # Test at 0°C: p_sat should be ~611 Pa (Table 17.2)
    prob = NonlinearProblem(ssys, Dict(ssys.T => 273.15))
    sol = solve(prob)
    @test sol[ssys.p_sat_water] ≈ 611.0 atol = 2.0  # ~611 Pa at 0°C

    # Test at 20°C: p_sat should be ~2338 Pa
    prob20 = NonlinearProblem(ssys, Dict(ssys.T => 293.15))
    sol20 = solve(prob20)
    @test sol20[ssys.p_sat_water] ≈ 2338.0 atol = 10.0  # ~2338 Pa at 20°C

    # Test at 100°C: p_sat should be ~101325 Pa (but polynomial valid only to 50°C)
    # Test at 30°C instead: should be ~4243 Pa
    prob30 = NonlinearProblem(ssys, Dict(ssys.T => 303.15))
    sol30 = solve(prob30)
    @test sol30[ssys.p_sat_water] ≈ 4243.0 atol = 20.0

    # Test ice at -10°C: p_sat_ice should be ~260 Pa
    prob_ice = NonlinearProblem(ssys, Dict(ssys.T => 263.15))
    sol_ice = solve(prob_ice)
    @test sol_ice[ssys.p_sat_ice] ≈ 260.0 atol = 5.0

    # Monotonicity: p_sat should increase with temperature
    @test sol20[ssys.p_sat_water] > sol[ssys.p_sat_water]
    @test sol30[ssys.p_sat_water] > sol20[ssys.p_sat_water]
end

@testitem "WaterProperties surface tension" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = WaterProperties()
    ssys = mtkcompile(sys)

    # Surface tension at 20°C: ~0.0728 N/m (Table 17.1)
    prob = NonlinearProblem(ssys, Dict(ssys.T => 293.15))
    sol = solve(prob)
    @test sol[ssys.σ_w0] ≈ 0.0730 atol = 0.002

    # Surface tension decreases with temperature
    prob0 = NonlinearProblem(ssys, Dict(ssys.T => 273.15))
    sol0 = solve(prob0)
    @test sol0[ssys.σ_w0] > sol[ssys.σ_w0]
end

@testitem "WaterProperties latent heat" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = WaterProperties()
    ssys = mtkcompile(sys)

    # Latent heat of vaporization at 20°C: ~2.45e6 J/kg (Table 17.1 value ~2.45 kJ/g)
    prob = NonlinearProblem(ssys, Dict(ssys.T => 293.15))
    sol = solve(prob)
    @test sol[ssys.ΔH_v] ≈ 2.45e6 atol = 0.05e6

    # Latent heat of melting at 0°C: ~333500 J/kg
    prob0 = NonlinearProblem(ssys, Dict(ssys.T => 273.15))
    sol0 = solve(prob0)
    @test sol0[ssys.ΔH_m] ≈ 333500.0 atol = 1000.0
end

@testitem "KelvinEffect structural" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = KelvinEffect()
    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 2
end

@testitem "KelvinEffect physical behavior" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = KelvinEffect()
    ssys = mtkcompile(sys)

    # For D_p = 0.1 µm, Kelvin ratio should be > 1 (enhanced vapor pressure)
    prob = NonlinearProblem(ssys, Dict(ssys.D_p => 1e-7, ssys.T => 293.15))
    sol = solve(prob)
    @test sol[ssys.kelvin_ratio] > 1.0
    @test sol[ssys.kelvin_ratio] < 1.1  # Should be small enhancement

    # For D_p = 1 µm, Kelvin effect should be negligible
    prob_large = NonlinearProblem(ssys, Dict(ssys.D_p => 1e-6, ssys.T => 293.15))
    sol_large = solve(prob_large)
    @test sol_large[ssys.kelvin_ratio] ≈ 1.0 atol = 0.003

    # Smaller droplets should have larger Kelvin effect
    prob_small = NonlinearProblem(ssys, Dict(ssys.D_p => 5e-8, ssys.T => 293.15))
    sol_small = solve(prob_small)
    @test sol_small[ssys.kelvin_ratio] > sol[ssys.kelvin_ratio]

    # Kelvin parameter A: should be ~2.1e-9 m at 293 K (Eq. 17.28)
    @test sol[ssys.A_kelvin] ≈ 2.1e-9 atol = 0.2e-9
end

@testitem "KohlerTheory structural" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = KohlerTheory()
    @test length(equations(sys)) == 9
    @test length(unknowns(sys)) == 9
end

@testitem "KohlerTheory critical supersaturation NaCl" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = KohlerTheory()
    ssys = mtkcompile(sys)

    # NaCl: M_s = 0.05844 kg/mol, ρ_s = 2165 kg/m³, ν = 2
    # For d_s = 0.1 µm, critical supersaturation should be ~0.1-0.3%
    prob = NonlinearProblem(ssys, Dict(
        ssys.d_s => 1e-7,    # 0.1 µm dry diameter
        ssys.M_s => 0.05844, # NaCl
        ssys.ρ_s => 2165.0,  # NaCl density
        ssys.ν_s => 2.0,     # NaCl dissociation
        ssys.T => 293.15,
        ssys.ε_m => 1.0,     # Fully soluble
    ))
    sol = solve(prob)

    # Critical supersaturation for 0.1 µm NaCl should be ~0.13%
    sc_pct = (sol[ssys.S_c] - 1) * 100
    @test sc_pct > 0.05  # At least 0.05%
    @test sc_pct < 0.5   # Less than 0.5%

    # Critical diameter should be larger than dry diameter
    @test sol[ssys.D_pc] > 1e-7

    # B parameter should be positive
    @test sol[ssys.B] > 0

    # Smaller dry particles should have higher critical supersaturation
    prob_small = NonlinearProblem(ssys, Dict(
        ssys.d_s => 5e-8,    # 0.05 µm dry diameter
        ssys.M_s => 0.05844,
        ssys.ρ_s => 2165.0,
        ssys.ν_s => 2.0,
        ssys.T => 293.15,
        ssys.ε_m => 1.0,
    ))
    sol_small = solve(prob_small)
    @test sol_small[ssys.S_c] > sol[ssys.S_c]  # Higher critical S for smaller particles
end

@testitem "KohlerTheory critical supersaturation (NH4)2SO4" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = KohlerTheory()
    ssys = mtkcompile(sys)

    # (NH4)2SO4: M_s = 0.13214 kg/mol, ρ_s = 1770 kg/m³, ν = 3
    prob = NonlinearProblem(ssys, Dict(
        ssys.d_s => 1e-7,     # 0.1 µm dry diameter
        ssys.M_s => 0.13214,  # (NH4)2SO4
        ssys.ρ_s => 1770.0,   # (NH4)2SO4 density
        ssys.ν_s => 3.0,      # (NH4)2SO4 dissociation
        ssys.T => 293.15,
        ssys.ε_m => 1.0,
    ))
    sol = solve(prob)

    # Critical supersaturation should be comparable to NaCl (slightly different)
    sc_pct = (sol[ssys.S_c] - 1) * 100
    @test sc_pct > 0.05
    @test sc_pct < 0.5
end

@testitem "KohlerTheory with insoluble material" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = KohlerTheory()
    ssys = mtkcompile(sys)

    # Fully soluble NaCl particle
    prob_full = NonlinearProblem(ssys, Dict(
        ssys.d_s => 1e-7, ssys.M_s => 0.05844, ssys.ρ_s => 2165.0,
        ssys.ν_s => 2.0, ssys.T => 293.15, ssys.ε_m => 1.0,
    ))
    sol_full = solve(prob_full)

    # Half soluble particle (ε_m = 0.5)
    prob_half = NonlinearProblem(ssys, Dict(
        ssys.d_s => 1e-7, ssys.M_s => 0.05844, ssys.ρ_s => 2165.0,
        ssys.ν_s => 2.0, ssys.T => 293.15, ssys.ε_m => 0.5,
    ))
    sol_half = solve(prob_half)

    # More insoluble material should increase critical supersaturation
    @test sol_half[ssys.S_c] > sol_full[ssys.S_c]
end

@testitem "DropletGrowth structural" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = DropletGrowth()
    @test length(equations(sys)) == 8
    @test length(unknowns(sys)) == 8
end

@testitem "DropletGrowth physical behavior" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = DropletGrowth()
    ssys = mtkcompile(sys)

    # In supersaturated environment (S_v > 1), droplet should grow
    prob = NonlinearProblem(ssys, Dict(
        ssys.S_v => 1.01, ssys.D_p => 1e-6, ssys.T => 293.15,
        ssys.n_s => 0.0, ssys.d_u => 0.0,
    ))
    sol = solve(prob)
    @test sol[ssys.dDp_dt] > 0  # Growing

    # In subsaturated environment (S_v < 1), droplet should shrink
    prob_sub = NonlinearProblem(ssys, Dict(
        ssys.S_v => 0.99, ssys.D_p => 1e-6, ssys.T => 293.15,
        ssys.n_s => 0.0, ssys.d_u => 0.0,
    ))
    sol_sub = solve(prob_sub)
    @test sol_sub[ssys.dDp_dt] < 0  # Shrinking

    # Diffusivity should be ~2.4e-5 m²/s at 293 K, 1 atm
    @test sol[ssys.D_v] ≈ 2.4e-5 atol = 0.3e-5

    # Thermal conductivity at 293 K: ~0.025 W/(m·K)
    @test sol[ssys.k_a] ≈ 0.025 atol = 0.003

    # Growth rate proportional to 1/D_p: smaller droplet should grow faster
    prob_small = NonlinearProblem(ssys, Dict(
        ssys.S_v => 1.01, ssys.D_p => 5e-7, ssys.T => 293.15,
        ssys.n_s => 0.0, ssys.d_u => 0.0,
    ))
    sol_small = solve(prob_small)
    @test sol_small[ssys.dDp_dt] > sol[ssys.dDp_dt]
end

@testitem "CloudDynamics structural" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = CloudDynamics()
    @test length(equations(sys)) == 5
    @test length(unknowns(sys)) == 5
end

@testitem "CloudDynamics dew point and LCL" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = CloudDynamics()
    ssys = mtkcompile(sys)

    # At 100% RH, dew point should equal temperature
    prob_100 = NonlinearProblem(ssys, Dict(ssys.T_0 => 293.15, ssys.RH => 1.0))
    sol_100 = solve(prob_100)
    @test sol_100[ssys.T_d] ≈ 293.15 atol = 0.1
    @test sol_100[ssys.h_LCL] ≈ 0.0 atol = 10.0

    # At lower RH, dew point should be below temperature
    prob_80 = NonlinearProblem(ssys, Dict(ssys.T_0 => 293.15, ssys.RH => 0.8))
    sol_80 = solve(prob_80)
    @test sol_80[ssys.T_d] < 293.15

    # LCL should be positive and increase with decreasing RH
    @test sol_80[ssys.h_LCL] > 0
    prob_60 = NonlinearProblem(ssys, Dict(ssys.T_0 => 293.15, ssys.RH => 0.6))
    sol_60 = solve(prob_60)
    @test sol_60[ssys.h_LCL] > sol_80[ssys.h_LCL]

    # Dry adiabatic lapse rate should be ~9.76 K/km = 0.00976 K/m
    @test sol_80[ssys.Γ_d] ≈ 0.00976 atol = 0.001
end

@testitem "IcePhysics structural" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = IcePhysics()
    @test length(equations(sys)) == 4
    @test length(unknowns(sys)) == 4
end

@testitem "IcePhysics physical behavior" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = IcePhysics()
    ssys = mtkcompile(sys)

    # Freezing point depression for pure water (a_w = 1): should be ~0
    prob = NonlinearProblem(ssys, Dict(ssys.T => 263.15, ssys.a_w => 1.0, ssys.molality => 0.0))
    sol = solve(prob)
    @test sol[ssys.ΔT_f] ≈ 0.0 atol = 0.01

    # For 1 mol/kg NaCl solution, ΔT_f ≈ 1.86 K (cryoscopic constant)
    prob_salt = NonlinearProblem(ssys, Dict(ssys.T => 263.15, ssys.a_w => 1.0, ssys.molality => 1.0))
    sol_salt = solve(prob_salt)
    @test sol_salt[ssys.ΔT_f] ≈ 1.86 atol = 0.1

    # Kelvin effect for ice: should be > 1 for small particles
    @test sol[ssys.kelvin_ice] > 1.0

    # Ice nuclei concentration increases with decreasing temperature
    prob_cold = NonlinearProblem(ssys, Dict(ssys.T => 253.15, ssys.a_w => 1.0))
    sol_cold = solve(prob_cold)
    prob_warm = NonlinearProblem(ssys, Dict(ssys.T => 263.15, ssys.a_w => 1.0))
    sol_warm = solve(prob_warm)
    @test sol_cold[ssys.IN] > sol_warm[ssys.IN]
end

@testitem "RainFormation structural" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = RainFormation()
    @test length(equations(sys)) == 4
    @test length(unknowns(sys)) == 4
end

@testitem "RainFormation physical behavior" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = RainFormation()
    ssys = mtkcompile(sys)

    prob = NonlinearProblem(ssys, Dict())
    sol = solve(prob)

    # Collection efficiency should be product of collision and coalescence
    @test sol[ssys.E_t] ≈ 0.5 * 1.0

    # Mass accretion rate should be positive
    @test sol[ssys.dm_dt] > 0

    # Marshall-Palmer distribution: n(D_p) should decrease with size
    prob_large = NonlinearProblem(ssys, Dict(ssys.D_p => 2e-3))
    sol_large = solve(prob_large)
    @test sol_large[ssys.n_MP] < sol[ssys.n_MP]

    # CDF should be between 0 and 1
    @test 0 ≤ sol[ssys.F_dist] ≤ 1
end

@testitem "AerosolScavenging structural" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = AerosolScavenging()
    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3
end

@testitem "AerosolScavenging physical behavior" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = AerosolScavenging()
    ssys = mtkcompile(sys)

    prob = NonlinearProblem(ssys, Dict())
    sol = solve(prob)

    # Mass scavenging ratio: (1e9 - 5e8)/1e9 = 0.5
    @test sol[ssys.F_mass] ≈ 0.5 atol = 0.01

    # Number scavenging ratio: (1e9 - 3e8)/1e9 = 0.7
    @test sol[ssys.F_number] ≈ 0.7 atol = 0.01

    # Scavenging coefficient: 1e8 * 1e-12 = 1e-4 s^-1
    @test sol[ssys.Λ] ≈ 1e-4 atol = 1e-6

    # Scavenging ratios should be between 0 and 1
    @test 0 ≤ sol[ssys.F_mass] ≤ 1
    @test 0 ≤ sol[ssys.F_number] ≤ 1
end

@testitem "CloudPhysics composite structural" tags=[:cloud_physics] begin
    using ModelingToolkit
    using Aerosol

    sys = CloudPhysics()
    @test length(equations(sys)) == 42  # 7+2+9+8+5+4+4+3

    # Verify all subsystems are present
    subsys_names = Symbol.(nameof.(ModelingToolkit.get_systems(sys)))
    for name in [:water, :kelvin, :kohler, :growth, :cloud, :ice, :rain, :scav]
        @test name in subsys_names
    end
end
