# Organic Atmospheric Aerosols

## Overview

This module implements key equations from Chapter 14 of Seinfeld and Pandis (2006) for
organic atmospheric aerosol modeling, including:

- The EC tracer method for estimating secondary organic carbon
- Gas-particle partitioning models for SOA formation (noninteracting, absorptive, two-product)
- Adsorption isotherms (Langmuir, BET, FHH)

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 14. John Wiley & Sons, Inc.

```@docs
ECTracerMethod
NoninteractingSOA
AbsorptivePartitioning
TwoProductSOA
LangmuirAdsorption
BETAdsorption
FHHAdsorption
```

## Implementation

### EC Tracer Method (Section 14.3.2)

The EC tracer method estimates secondary organic carbon from measurements of total OC and EC.

```@example organic_aerosol
using Aerosol
using ModelingToolkit
using DataFrames, Symbolics, DynamicQuantities

sys = ECTracerMethod()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example organic_aerosol
equations(sys)
```

### Noninteracting SOA (Section 14.5.2)

Gas-particle partitioning for a compound that condenses as a pure phase.

```@example organic_aerosol
sys_ni = NoninteractingSOA()
vars = unknowns(sys_ni)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example organic_aerosol
equations(sys_ni)
```

### Absorptive Partitioning (Section 14.5.2)

SOA formation by absorption into a preexisting organic aerosol phase.

```@example organic_aerosol
sys_ap = AbsorptivePartitioning()
vars = unknowns(sys_ap)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example organic_aerosol
equations(sys_ap)
```

### Two-Product SOA Model (Section 14.5.2)

The Odum et al. (1996) two-product model for SOA formation.

```@example organic_aerosol
sys_tp = TwoProductSOA()
vars = unknowns(sys_tp)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example organic_aerosol
equations(sys_tp)
```

### Adsorption Isotherms (Section 14.5.3)

```@example organic_aerosol
sys_lang = LangmuirAdsorption()
equations(sys_lang)
```

```@example organic_aerosol
sys_bet = BETAdsorption()
equations(sys_bet)
```

```@example organic_aerosol
sys_fhh = FHHAdsorption()
equations(sys_fhh)
```

## Analysis

### Noninteracting SOA: Aerosol Yield vs Reacted ROG

This figure demonstrates the threshold behavior of noninteracting SOA formation.
Below the threshold ΔROG, all product remains in the gas phase. Above it,
the aerosol yield increases with increasing reacted ROG.

```@example organic_aerosol
using NonlinearSolve, Plots

sys_ni = NoninteractingSOA()
sys_ni_nns = ModelingToolkit.toggle_namespacing(sys_ni, false)
ssys_ni = mtkcompile(sys_ni; inputs=[sys_ni_nns.ΔROG])

# Compute threshold from default parameters
R = 8.314; T = 298.0; p_i = 1.01325e-5; M_i = 0.180; M_ROG = 0.150; a_i = 0.05
threshold = p_i * M_ROG / (a_i * R * T)

# Sweep ΔROG from 0 to 20× threshold
ΔROG_vals = range(0, 20 * threshold, length=200)
Y_vals = Float64[]
for drog in ΔROG_vals
    prob = NonlinearProblem(ssys_ni, Dict(
        ssys_ni.ΔROG => drog,
        ssys_ni.c_eq => 1e-10, ssys_ni.c_total => 1e-10,
        ssys_ni.c_aer => 1e-10, ssys_ni.c_gas => 1e-10,
        ssys_ni.ΔROG_threshold => 1e-10,
        ssys_ni.X_p => 0.5, ssys_ni.Y => 0.01))
    sol = solve(prob)
    push!(Y_vals, sol[ssys_ni.Y])
end

plot(ΔROG_vals ./ threshold, Y_vals,
    xlabel="ΔROG / ΔROG*", ylabel="Aerosol Mass Yield Y",
    title="Noninteracting SOA: Yield vs Reacted ROG",
    label="Y (Eq. 14.16)", legend=:bottomright, lw=2)
vline!([1.0], label="Threshold ΔROG*", ls=:dash, color=:red)
```

### Absorptive Partitioning: Effect of Preexisting OA

When preexisting organic aerosol is present, there is no threshold for SOA formation.
The aerosol mass fraction `X_p` depends on the ratio of absorbing medium concentration
to vapor pressure.

```@example organic_aerosol
sys_ap = AbsorptivePartitioning()
sys_ap_nns = ModelingToolkit.toggle_namespacing(sys_ap, false)
ssys_ap = mtkcompile(sys_ap; inputs=[sys_ap_nns.ΔROG])

# Sweep m_0 (preexisting OA) from 0.1 to 100 μg/m³
m0_vals = 10 .^ range(-10, -7, length=100)  # kg/m³
Xp_vals = Float64[]
for m0 in m0_vals
    prob = NonlinearProblem(ssys_ap, Dict(
        ssys_ap.ΔROG => 50e-9,
        ssys_ap.c_aer => 1e-10, ssys_ap.X_p => 0.5, ssys_ap.Y => 0.01,
        ssys_ap.m_0 => m0))
    sol = solve(prob)
    push!(Xp_vals, sol[ssys_ap.X_p])
end

plot(m0_vals * 1e9, Xp_vals,
    xlabel="Preexisting OA m₀ (μg/m³)", ylabel="Mass Fraction in Aerosol X_p",
    title="Absorptive Partitioning: X_p vs m₀ (Eq. 14.25)",
    label="X_p", xscale=:log10, lw=2, legend=:bottomright)
```

### Langmuir Adsorption Isotherm

The Langmuir isotherm shows saturation behavior at high pressures.

```@example organic_aerosol
sys_lang = LangmuirAdsorption()
sys_lang_nns = ModelingToolkit.toggle_namespacing(sys_lang, false)
ssys_lang = mtkcompile(sys_lang; inputs=[sys_lang_nns.p])

p_vals = 10 .^ range(-3, 3, length=200)
V_vals = Float64[]
for pv in p_vals
    prob = NonlinearProblem(ssys_lang, Dict(ssys_lang.p => pv, ssys_lang.V => 0.5))
    sol = solve(prob)
    push!(V_vals, sol[ssys_lang.V])
end

plot(p_vals, V_vals,
    xlabel="Pressure p (Pa)", ylabel="Adsorbed Volume V (m³)",
    title="Langmuir Adsorption Isotherm (Eq. 14.44)",
    label="V/V_m", xscale=:log10, lw=2, legend=:bottomright)
hline!([1.0], label="V_m (monolayer)", ls=:dash, color=:red)
```

### Comparison of BET and Langmuir Isotherms

The BET isotherm extends the Langmuir model to multilayer adsorption.
At low saturation ratios, both isotherms agree; at higher values,
BET predicts additional adsorption layers.

```@example organic_aerosol
sys_bet = BETAdsorption()
sys_bet_nns = ModelingToolkit.toggle_namespacing(sys_bet, false)
ssys_bet = mtkcompile(sys_bet; inputs=[sys_bet_nns.S])

S_vals = range(0.01, 0.95, length=200)
V_bet = Float64[]
V_lang = Float64[]
c_BET = 10.0

for sv in S_vals
    # BET
    prob_bet = NonlinearProblem(ssys_bet, Dict(ssys_bet.S => sv, ssys_bet.V => 1.0))
    sol_bet = solve(prob_bet)
    push!(V_bet, sol_bet[ssys_bet.V])

    # Equivalent Langmuir: b*p = c_BET*S
    prob_lang = NonlinearProblem(ssys_lang, Dict(ssys_lang.p => c_BET * sv, ssys_lang.V => 0.5))
    sol_lang = solve(prob_lang)
    push!(V_lang, sol_lang[ssys_lang.V])
end

plot(S_vals, V_bet, label="BET (c=10)", lw=2,
    xlabel="Saturation Ratio S", ylabel="V / V_m",
    title="BET vs Langmuir Isotherms")
plot!(S_vals, V_lang, label="Langmuir (b·p = c·S)", lw=2, ls=:dash)
hline!([1.0], label="Monolayer", ls=:dot, color=:gray)
```
