# Cloud Physics

## Overview

This module implements the cloud microphysics equations from Chapter 17 of Seinfeld & Pandis (2006), covering thermodynamic properties of water, the Kelvin (curvature) effect, Köhler theory for cloud droplet activation, diffusional growth of cloud droplets, adiabatic cooling and cloud formation dynamics, ice nucleation, rain formation, and aerosol scavenging by cloud droplets.

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 2nd Edition, John Wiley & Sons, Inc. Chapter 17: Cloud Physics.

```@docs
WaterProperties
KelvinEffect
KohlerTheory
DropletGrowth
CloudDynamics
IcePhysics
RainFormation
AerosolScavenging
CloudPhysics
```

## Implementation

The cloud physics system is composed of eight subsystems, each implemented as a ModelingToolkit `@component` function. They can be used individually or composed together via the `CloudPhysics` parent system.

### State Variables

```@example cloud_physics
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Aerosol

sys = CloudPhysics()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example cloud_physics
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example cloud_physics
eqs = equations(sys)
```

## Analysis

### Figure 17.1: Saturation Vapor Pressure vs Temperature

The saturation vapor pressure of water increases approximately exponentially with temperature (Table 17.2 polynomial).

```@example cloud_physics
using Plots, NonlinearSolve

wp = WaterProperties()
wp_sys = mtkcompile(wp)

T_range = 263.15:1.0:313.15  # -10°C to 40°C
p_sat_vals = Float64[]

for T_val in T_range
    prob = NonlinearProblem(wp_sys, Dict(wp_sys.T => T_val))
    sol = solve(prob)
    push!(p_sat_vals, sol[wp_sys.p_sat_water])
end

plot(T_range .- 273.15, p_sat_vals ./ 100,
    xlabel="Temperature (°C)", ylabel="Saturation Vapor Pressure (hPa)",
    label="Water (Table 17.2 polynomial)",
    linewidth=2, legend=:topleft,
    title="Fig. 17.1: Saturation Vapor Pressure")
```

### Figure 17.2: Kelvin Effect

The vapor pressure over a curved droplet surface is enhanced relative to a flat surface, with the enhancement increasing as droplet diameter decreases (Eq. 17.9).

```@example cloud_physics
ke = KelvinEffect()
ke_sys = mtkcompile(ke)

D_p_range = 10.0 .^ range(log10(1e-8), log10(1e-5), length=100)  # 10 nm to 10 µm
kelvin_vals = Float64[]

for D_val in D_p_range
    prob = NonlinearProblem(ke_sys, Dict(ke_sys.D_p => D_val, ke_sys.T => 293.15))
    sol = solve(prob)
    push!(kelvin_vals, sol[ke_sys.kelvin_ratio])
end

plot(D_p_range * 1e6, kelvin_vals,
    xlabel="Droplet Diameter (µm)", ylabel="p_w(D_p) / p°",
    xscale=:log10, label="T = 293.15 K",
    linewidth=2, title="Fig. 17.2: Kelvin Effect",
    legend=:topright)
hline!([1.0], linestyle=:dash, color=:gray, label="Flat surface")
```

### Figure 17.5: Köhler Curves

Köhler curves show the equilibrium saturation ratio as a function of wet droplet diameter for different dry particle sizes and solute compositions. The curves exhibit a maximum (critical supersaturation) beyond which droplets activate and grow without bound.

```@example cloud_physics
kt = KohlerTheory()
kt_sys = mtkcompile(kt)

D_wet_range = 10.0 .^ range(log10(1e-8), log10(5e-5), length=200)

function kohler_curve(kt_sys, d_dry, M_s, rho_s, nu_s)
    S_vals = Float64[]
    for D_w in D_wet_range
        prob = NonlinearProblem(kt_sys, Dict(
            kt_sys.D_p => D_w, kt_sys.d_s => d_dry,
            kt_sys.M_s => M_s, kt_sys.ρ_s => rho_s,
            kt_sys.ν_s => nu_s, kt_sys.T => 293.15, kt_sys.ε_m => 1.0
        ))
        sol = solve(prob)
        push!(S_vals, sol[kt_sys.S])
    end
    return S_vals
end

# NaCl: M_s = 0.05844, rho = 2165, nu = 2
# (NH4)2SO4: M_s = 0.13214, rho = 1770, nu = 3
plot(title="Fig. 17.5: Köhler Curves", xlabel="Wet Diameter (µm)",
    ylabel="Saturation Ratio S", legend=:topright, ylims=(0.990, 1.015))

for (d_dry, style) in [(0.05e-6, :solid), (0.1e-6, :dash), (0.5e-6, :dot)]
    d_label = round(d_dry * 1e6, digits=2)
    S_nacl = kohler_curve(kt_sys, d_dry, 0.05844, 2165.0, 2.0)
    plot!(D_wet_range * 1e6, S_nacl, label="NaCl d_s=$(d_label) µm",
        linestyle=style, linewidth=2, color=:blue, xscale=:log10)
    S_as = kohler_curve(kt_sys, d_dry, 0.13214, 1770.0, 3.0)
    plot!(D_wet_range * 1e6, S_as, label="(NH₄)₂SO₄ d_s=$(d_label) µm",
        linestyle=style, linewidth=2, color=:red, xscale=:log10)
end
hline!([1.0], linestyle=:dash, color=:gray, label="S = 1")
plot!()
```

### Figure 17.6: Critical Supersaturation vs Dry Diameter

Critical supersaturation decreases with increasing dry particle diameter, following a power-law relationship (Eq. 17.34).

```@example cloud_physics
d_dry_range = 10.0 .^ range(log10(1e-8), log10(1e-6), length=50)

function critical_ss(kt_sys, d_dry_range, M_s, rho_s, nu_s)
    sc_vals = Float64[]
    for d_d in d_dry_range
        prob = NonlinearProblem(kt_sys, Dict(
            kt_sys.d_s => d_d, kt_sys.M_s => M_s, kt_sys.ρ_s => rho_s,
            kt_sys.ν_s => nu_s, kt_sys.T => 293.15, kt_sys.ε_m => 1.0
        ))
        sol = solve(prob)
        push!(sc_vals, (sol[kt_sys.S_c] - 1) * 100)  # Convert to percent
    end
    return sc_vals
end

sc_nacl = critical_ss(kt_sys, d_dry_range, 0.05844, 2165.0, 2.0)
sc_as = critical_ss(kt_sys, d_dry_range, 0.13214, 1770.0, 3.0)

plot(d_dry_range * 1e6, sc_nacl, label="NaCl", linewidth=2, color=:blue,
    xlabel="Dry Diameter (µm)", ylabel="Critical Supersaturation (%)",
    xscale=:log10, yscale=:log10, title="Fig. 17.6: Critical Supersaturation",
    legend=:topright)
plot!(d_dry_range * 1e6, sc_as, label="(NH₄)₂SO₄", linewidth=2, color=:red)
```

### Figure 17.9: Effect of Insoluble Material on Critical Supersaturation

Particles with less soluble material require higher supersaturation to activate (Eqs. 17.35-17.40).

```@example cloud_physics
epsilon_vals = [1.0, 0.5, 0.1, 0.01]

plot(title="Fig. 17.9: Insoluble Material Effect",
    xlabel="Dry Diameter (µm)", ylabel="Critical Supersaturation (%)",
    xscale=:log10, yscale=:log10, legend=:topright)

for eps_m in epsilon_vals
    sc_vals = Float64[]
    for d_d in d_dry_range
        prob = NonlinearProblem(kt_sys, Dict(
            kt_sys.d_s => d_d, kt_sys.M_s => 0.13214, kt_sys.ρ_s => 1770.0,
            kt_sys.ν_s => 3.0, kt_sys.T => 293.15, kt_sys.ε_m => eps_m,
            kt_sys.ρ_u => 2000.0
        ))
        sol = solve(prob)
        push!(sc_vals, (sol[kt_sys.S_c] - 1) * 100)
    end
    plot!(d_dry_range * 1e6, sc_vals, label="ε_m = $eps_m", linewidth=2)
end
plot!()
```

### Figure 17.10: Dew Point Temperature

Dew point temperature as a function of ambient temperature for different relative humidities (Eq. 17.43).

```@example cloud_physics
cd = CloudDynamics()
cd_sys = mtkcompile(cd)

T0_range = 273.15:2.0:313.15  # 0°C to 40°C
RH_vals = [0.2, 0.4, 0.6, 0.8, 1.0]

plot(title="Fig. 17.10: Dew Point Temperature",
    xlabel="Temperature (°C)", ylabel="Dew Point Temperature (°C)",
    legend=:bottomright)

for rh in RH_vals
    Td_vals = Float64[]
    for T0 in T0_range
        prob = NonlinearProblem(cd_sys, Dict(cd_sys.T_0 => T0, cd_sys.RH => rh))
        sol = solve(prob)
        push!(Td_vals, sol[cd_sys.T_d] - 273.15)
    end
    plot!(T0_range .- 273.15, Td_vals, label="RH = $rh", linewidth=2)
end
plot!()
```

### Figure 17.11: Lifting Condensation Level

Height of the lifting condensation level (LCL) as a function of surface temperature and relative humidity (Eq. 17.49).

```@example cloud_physics
plot(title="Fig. 17.11: Lifting Condensation Level",
    xlabel="Surface Temperature (°C)", ylabel="LCL Height (m)",
    legend=:topleft)

for rh in [0.2, 0.4, 0.6, 0.8]
    h_vals = Float64[]
    for T0 in T0_range
        prob = NonlinearProblem(cd_sys, Dict(cd_sys.T_0 => T0, cd_sys.RH => rh))
        sol = solve(prob)
        push!(h_vals, sol[cd_sys.h_LCL])
    end
    plot!(T0_range .- 273.15, h_vals, label="RH = $rh", linewidth=2)
end
plot!()
```

### Droplet Growth Rate

Growth rate of cloud droplets as a function of diameter at 1% supersaturation (Eq. 17.70), demonstrating that smaller droplets grow faster (dD_p/dt ∝ 1/D_p).

```@example cloud_physics
dg = DropletGrowth()
dg_sys = mtkcompile(dg)

D_range = 10.0 .^ range(log10(1e-6), log10(50e-6), length=50)
growth_vals = Float64[]

for D_val in D_range
    prob = NonlinearProblem(dg_sys, Dict(
        dg_sys.D_p => D_val, dg_sys.S_v => 1.01, dg_sys.T => 293.15,
        dg_sys.n_s => 0.0, dg_sys.d_u => 0.0
    ))
    sol = solve(prob)
    push!(growth_vals, sol[dg_sys.dDp_dt])
end

plot(D_range * 1e6, growth_vals * 1e6,
    xlabel="Droplet Diameter (µm)", ylabel="Growth Rate dD_p/dt (µm/s)",
    linewidth=2, label="S = 1.01, T = 293 K",
    title="Droplet Growth Rate (Eq. 17.70)",
    legend=:topright)
```

### Ice Nuclei Concentration

Ice nuclei concentration increases exponentially with decreasing temperature (Eq. 17.103).

```@example cloud_physics
ice = IcePhysics()
ice_sys = mtkcompile(ice)

T_ice_range = 233.15:1.0:273.15  # -40°C to 0°C
IN_vals = Float64[]

for T_val in T_ice_range
    prob = NonlinearProblem(ice_sys, Dict(ice_sys.T => T_val, ice_sys.a_w => 1.0))
    sol = solve(prob)
    push!(IN_vals, sol[ice_sys.IN] / 1000.0)  # Convert from /m³ to /L
end

plot(T_ice_range .- 273.15, IN_vals,
    xlabel="Temperature (°C)", ylabel="Ice Nuclei (L⁻¹)",
    yscale=:log10, linewidth=2, label="Eq. 17.103",
    title="Ice Nuclei Concentration",
    legend=:topright)
```
