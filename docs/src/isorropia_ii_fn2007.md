# ISORROPIA II (Fountoukis & Nenes, 2007)

## Overview

ISORROPIA II is an inorganic aerosol thermodynamic equilibrium model for the
K⁺–Ca²⁺–Mg²⁺–NH₄⁺–Na⁺–SO₄²⁻–NO₃⁻–Cl⁻–H₂O system. It computes the equilibrium
partitioning of semi-volatile inorganic species between gas, aqueous, and solid phases.

The implementation supports two solution modes via the `stable` parameter:
- **Metastable** (`stable=0`, default): No solid precipitation — the aerosol is assumed to be a liquid solution at all relative humidities.
- **Stable** (`stable=1`): Solid precipitation is allowed — salts crystallize when the solution becomes supersaturated. The mutual deliquescence relative humidity (MDRH) from Table 5 is enforced: below the MDRH for the given salt mixture, the aerosol is completely dry (crystalline). Temperature-dependent DRH values follow Eq. 17 and Table 4.

**Reference**: Fountoukis, C. and Nenes, A.: ISORROPIA II: a computationally
efficient thermodynamic equilibrium model for
K⁺–Ca²⁺–Mg²⁺–NH₄⁺–Na⁺–SO₄²⁻–NO₃⁻–Cl⁻–H₂O aerosols,
Atmos. Chem. Phys., 7, 4639–4659, 2007.

```@docs
IsorropiaEquilibrium
```

## Implementation

The model solves a system of coupled algebraic equations that depends on the solution mode:

**Metastable mode** (`stable=0`, default): 21 equations total
- 8 mass balance equations (no solid terms)
- 1 electroneutrality (charge balance) equation
- 5 aqueous equilibrium expressions (HSO₄⁻, NH₃, HNO₃, HCl, H₂O)
- 1 ZSR water content equation (Eq. 16)
- 1 ionic strength definition
- 5 activity coefficient equations (Kusik-Meissner, Eqs. 9-13, with temperature correction Eq. 14)

**Stable mode** (`stable=1`): Additional equations for solid precipitation
- All metastable equations plus 19+ solid equilibrium equations using complementarity constraints

The solid precipitation uses a smooth complementarity formulation: for each solid salt,
either its amount is zero (undersaturated) or the ion activity product equals the
solubility product (at dissolution equilibrium). In metastable mode (`stable=0`),
all solid equations reduce to s = 0.

### State Variables

```@example iso2
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Aerosol

sys = IsorropiaEquilibrium()
compiled = mtkcompile(sys)

vars = unknowns(compiled)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example iso2
params = parameters(compiled)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Default => [ModelingToolkit.getdefault(p) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example iso2
equations(sys)
```

## Analysis

### Equilibrium Constants vs Temperature (cf. Table 2)

The equilibrium constants follow the van't Hoff temperature dependence (Eq. 5):

```math
K(T) = K_0 \exp\left[A\left(\frac{T_0}{T} - 1\right) + B\left(1 + \ln\frac{T_0}{T} - \frac{T_0}{T}\right)\right]
```

```@example iso2
using Plots

T_range = collect(260.0:1.0:320.0)

K1_vals = [Aerosol._iso2_eq_const(1.015e-2, 8.85, 25.14, T) for T in T_range]
K4_vals = [Aerosol._iso2_eq_const(2.511e6, 29.17, 16.83, T) for T in T_range]
K3_vals = [Aerosol._iso2_eq_const(1.971e6, 30.20, 19.91, T) for T in T_range]
Kw_vals = [Aerosol._iso2_eq_const(1.010e-14, -22.52, 26.92, T) for T in T_range]

p = plot(T_range, K1_vals, label="K1 (HSO4 dissoc.)", yscale=:log10,
    xlabel="Temperature (K)", ylabel="Equilibrium Constant",
    title="Temperature Dependence of Equilibrium Constants (Eq. 5)")
plot!(p, T_range, K4_vals, label="K4 (HNO3)")
plot!(p, T_range, K3_vals, label="K3 (HCl)")
plot!(p, T_range, Kw_vals, label="Kw (H2O)")
p
```

### Activity Coefficients vs Ionic Strength (cf. Eqs. 9-13)

Kusik-Meissner binary mean activity coefficients at 298.15 K:

```@example iso2
I_range = collect(0.01:0.1:30.0)

gamma_NaCl = [Aerosol._iso2_km_gamma(2.23, 1, I) for I in I_range]
gamma_NH4Cl = [Aerosol._iso2_km_gamma(0.82, 1, I) for I in I_range]
gamma_H2SO4 = [Aerosol._iso2_km_gamma(-0.10, 2, I) for I in I_range]
gamma_HNO3 = [Aerosol._iso2_km_gamma(2.60, 1, I) for I in I_range]
gamma_HCl = [Aerosol._iso2_km_gamma(6.00, 1, I) for I in I_range]

p = plot(I_range, gamma_NaCl, label="NaCl (q=2.23)", yscale=:log10,
    xlabel="Ionic Strength (mol/kg)", ylabel="Mean Activity Coefficient",
    title="Kusik-Meissner Activity Coefficients (Eqs. 9-13)")
plot!(p, I_range, gamma_NH4Cl, label="NH4Cl (q=0.82)")
plot!(p, I_range, gamma_H2SO4, label="H2SO4 (q=-0.10)")
plot!(p, I_range, gamma_HNO3, label="HNO3 (q=2.60)")
plot!(p, I_range, gamma_HCl, label="HCl (q=6.00)")
p
```

### ZSR Water Content vs RH (cf. Eq. 16)

Aerosol water content for pure salt aerosols computed using the ZSR method:

```@example iso2
RH_range = collect(0.10:0.01:0.98)

# (NH4)2SO4: c_NH4 = 2e-7 mol/m³, c_SO4 = 1e-7 mol/m³
W_as = [Aerosol._iso2_zsr_water(rh, 0.0, 2e-7, 1e-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        for rh in RH_range]

# NaCl: c_Na = 1e-7 mol/m³, c_Cl = 1e-7 mol/m³
W_nacl = [Aerosol._iso2_zsr_water(rh, 1e-7, 0.0, 0.0, 0.0, 0.0, 1e-7, 0.0, 0.0, 0.0)
          for rh in RH_range]

p = plot(RH_range, W_as, label="(NH4)2SO4",
    xlabel="Relative Humidity", ylabel="Water Content (kg/m3)",
    title="ZSR Aerosol Water Content (Eq. 16)")
plot!(p, RH_range, W_nacl, label="NaCl")
p
```

### DRH vs Temperature (cf. Figure 2, Eq. 17, Table 4)

Deliquescence relative humidity (DRH) for individual salts as a function of temperature.
DRH(T) = DRH(T₀) × exp(coeff × (1/T - 1/T₀)) where T₀ = 298.15 K. Salts with
large positive coefficients (e.g., NH₄NO₃, Ca(NO₃)₂) show strong temperature
dependence, while salts with near-zero coefficients (e.g., CaSO₄, KHSO₄, KNO₃)
are essentially temperature-independent.

```@example iso2
T_drh = collect(250.0:1.0:320.0)

# Select representative salts spanning different behaviors
salts = [
    (:NH42SO4, "(NH₄)₂SO₄"), (:NH4NO3, "NH₄NO₃"), (:NaCl, "NaCl"),
    (:Na2SO4, "Na₂SO₄"), (:NH4HSO4, "NH₄HSO₄"), (:CaNO32, "Ca(NO₃)₂"),
    (:K2SO4, "K₂SO₄"), (:MgSO4, "MgSO₄"), (:KCl, "KCl"),
]

p = plot(xlabel="Temperature (K)", ylabel="DRH",
    title="Deliquescence RH vs Temperature (Eq. 17, Table 4)",
    legend=:outerright, size=(800, 450))
for (sym, label) in salts
    drh0, coeff = Aerosol._ISO2_DRH[sym]
    drh_vals = [Aerosol._iso2_drh_T(drh0, coeff, T) for T in T_drh]
    plot!(p, T_drh, drh_vals, label=label)
end
p
```

### MDRH: Mutual Deliquescence (Table 5)

The mutual deliquescence relative humidity (MDRH) is the minimum RH at which any
aqueous phase can exist for a mixture of salts. Below MDRH, the aerosol is completely
crystalline. The MDRH depends on the salt mixture composition (sulfate ratio R₁ and
presence of crustal/sodium species). In stable mode, the water content is smoothly
driven to zero below MDRH using a logistic transition function.

| Composition Regime | MDRH |
|---|---|
| Crustal species present | 0.200 |
| Na present, no crustal (sulfate-poor) | 0.363 |
| NH₄-only (sulfate-poor) | 0.540 |
| NH₄-only (sulfate-rich) | 0.400 |
| Na present (sulfate-rich) | 0.363 |
| Crustal species present (sulfate-rich) | 0.200 |

### Equilibrium Partitioning: Gas-Aerosol Split

The fraction of semi-volatile species in the aerosol phase as a function of RH,
for a sulfate-poor continental aerosol (SO₄=1e-7, NH₃=3e-7, NO₃=1e-7 mol/m³):

```@example iso2
using NonlinearSolve, SciMLBase

RH_vals = collect(0.70:0.05:0.95)

NH4_frac = Float64[]
NO3_frac = Float64[]

for rh in RH_vals
    op = Pair{Symbolics.Num, Float64}[
        compiled.RH => rh,
        compiled.W_SO4_total => 1e-7,
        compiled.W_NH3_total => 3e-7,
        compiled.W_NO3_total => 1e-7,
        compiled.c_SO4 => 9.5e-8,
        compiled.c_NO3 => 5e-8,
        compiled.c_Cl => 1e-20,
        compiled.c_H => 1e-11,
        compiled.c_OH => 1e-14,
        compiled.I_s => 10.0,
    ]
    # Fill in default guesses for unknowns not explicitly set
    op_keys = Set(first.(op))
    for u in unknowns(compiled)
        if !(u in op_keys)
            g_val = ModelingToolkit.getguess(u)
            if g_val !== nothing
                push!(op, u => Float64(g_val))
            end
        end
    end
    prob = NonlinearProblem(compiled, op)
    sol = solve(prob, RobustMultiNewton(); maxiters=10000)
    if sol.retcode == SciMLBase.ReturnCode.Success
        push!(NH4_frac, sol[compiled.c_NH4] / 3e-7)
        push!(NO3_frac, sol[compiled.c_NO3] / 1e-7)
    else
        push!(NH4_frac, NaN)
        push!(NO3_frac, NaN)
    end
end

p = plot(RH_vals, NH4_frac, label="NH4+ / (NH4+ + NH3)",
    xlabel="Relative Humidity", ylabel="Aerosol Fraction",
    title="Gas-Aerosol Partitioning vs RH (Sulfate-Poor)")
plot!(p, RH_vals, NO3_frac, label="NO3- / (NO3- + HNO3)")
p
```

### Figures 6–9: Model Intercomparison Cases (Table 8)

The following figures reproduce Figures 6–9 from Fountoukis and Nenes (2007), showing
equilibrium aerosol composition vs. relative humidity for four representative cases from
Table 8. Input concentrations are converted from µg/m³ to mol/m³. The metastable
solution is shown (set `stable=1` to include solid precipitation).

```@example iso2
# Helper: sweep RH from high to low using continuation for solver convergence
function solve_iso2_sweep(compiled, Na, SO4, NH3, NO3, Cl, Ca, K, Mg;
                          RH_range=0.95:-0.01:0.30)
    RH_ok, W_vals = Float64[], Float64[]
    c = Dict(s => Float64[] for s in [:Na,:NH4,:SO4,:HSO4,:NO3,:Cl,:Ca,:K,:Mg,:H,:OH])
    g = Dict(s => Float64[] for s in [:NH3,:HNO3,:HCl])
    prev = nothing

    for rh in RH_range
        params = [compiled.RH => rh, compiled.W_Na_total => Na,
                  compiled.W_SO4_total => SO4, compiled.W_NH3_total => NH3,
                  compiled.W_NO3_total => NO3, compiled.W_Cl_total => Cl,
                  compiled.W_Ca_total => Ca, compiled.W_K_total => K,
                  compiled.W_Mg_total => Mg]
        if prev !== nothing
            guesses = [compiled.c_SO4 => prev[:SO4], compiled.c_HSO4 => prev[:HSO4],
                       compiled.c_NH4 => prev[:NH4], compiled.c_NO3 => prev[:NO3],
                       compiled.c_Cl => prev[:Cl], compiled.c_H => prev[:H],
                       compiled.c_OH => prev[:OH], compiled.I_s => prev[:I],
                       compiled.W_w => prev[:W]]
        else
            # Sulfate ratio R₁ determines aerosol regime
            R1 = SO4 > 0 ? (NH3 + Na + 2*Ca + K + 2*Mg) / SO4 : 100.0
            if R1 < 2  # Sulfate-rich: acidic aerosol
                guesses = [compiled.c_SO4 => max(0.5*SO4, 1e-20),
                           compiled.c_HSO4 => max(0.5*SO4, 1e-20),
                           compiled.c_NH4 => max(0.99*NH3, 1e-20),
                           compiled.c_NO3 => 1e-10, compiled.c_Cl => 1e-20,
                           compiled.c_H => 1e-7, compiled.c_OH => 1e-18,
                           compiled.I_s => 15.0, compiled.W_w => 1e-7,
                           compiled.g_NH3 => 1e-10, compiled.g_HNO3 => max(0.99*NO3, 1e-20)]
            else  # Sulfate-poor: more neutral aerosol
                guesses = [compiled.c_SO4 => max(0.85*SO4, 1e-20),
                           compiled.c_HSO4 => max(0.15*SO4, 1e-20),
                           compiled.c_NH4 => max(min(NH3, 2*SO4)*0.8, 1e-20),
                           compiled.c_NO3 => max(0.3*NO3, 1e-20),
                           compiled.c_Cl => max(0.3*Cl, 1e-20),
                           compiled.c_H => 1e-10, compiled.c_OH => 1e-13,
                           compiled.I_s => 10.0]
            end
        end

        # Fill in default guesses for any unknowns not explicitly set
        op = vcat(params, guesses)
        op_keys = Set(first.(op))
        for u in unknowns(compiled)
            if !(u in op_keys)
                g_val = ModelingToolkit.getguess(u)
                if g_val !== nothing
                    push!(op, u => Float64(g_val))
                end
            end
        end
        prob = NonlinearProblem(compiled, op)
        sol = solve(prob, RobustMultiNewton(); maxiters=10000)

        if sol.retcode == SciMLBase.ReturnCode.Success
            push!(RH_ok, rh); push!(W_vals, sol[compiled.W_w])
            for s in [:Na,:NH4,:SO4,:HSO4,:NO3,:Cl,:Ca,:K,:Mg,:H,:OH]
                push!(c[s], sol[getproperty(compiled, Symbol(:c_, s))])
            end
            push!(g[:NH3], sol[compiled.g_NH3])
            push!(g[:HNO3], sol[compiled.g_HNO3])
            push!(g[:HCl], sol[compiled.g_HCl])
            prev = Dict(:SO4=>sol[compiled.c_SO4], :HSO4=>sol[compiled.c_HSO4],
                        :NH4=>sol[compiled.c_NH4], :NO3=>sol[compiled.c_NO3],
                        :Cl=>sol[compiled.c_Cl], :H=>sol[compiled.c_H],
                        :OH=>sol[compiled.c_OH], :I=>sol[compiled.I_s],
                        :W=>sol[compiled.W_w])
        end
    end
    reverse!(RH_ok); reverse!(W_vals)
    for v in values(c); reverse!(v); end
    for v in values(g); reverse!(v); end
    return RH_ok, W_vals, c, g
end
nothing # hide
```

### Figure 6: Urban Aerosol — Sulfate Rich (Case 3, Table 8)

Aerosol composition vs. RH for an urban sulfate-rich case (1 < R₁ < 2):
Na=0, H₂SO₄=15.0, NH₃=2.0, HNO₃=10.0, HCl=0, Ca=0.9, K=1.0, Mg=0 µg/m³.
Panels show (a) aerosol water, (b) K⁺, (c) NH₄⁺, and (d) NO₃⁻, all in µg/m³.

```@example iso2
# Case 3 (Urban, sulfate rich) from Table 8
RH6, W6, c6, g6 = solve_iso2_sweep(compiled,
    0.0, 15.0e-6/98.08, 2.0e-6/17.03, 10.0e-6/63.01, 0.0,
    0.9e-6/40.08, 1.0e-6/39.10, 0.0)

p6a = plot(RH6, W6 .* 1e9, label="Metastable", xlabel="RH", ylabel="µg/m³",
    title="(a) Aerosol Water", legend=:topleft)
p6b = plot(RH6, c6[:K] .* 39.10e6, label="K⁺", xlabel="RH", ylabel="µg/m³",
    title="(b) Potassium")
p6c = plot(RH6, c6[:NH4] .* 18.04e6, label="NH₄⁺", xlabel="RH", ylabel="µg/m³",
    title="(c) Ammonium")
p6d = plot(RH6, c6[:NO3] .* 62.00e6, label="NO₃⁻", xlabel="RH", ylabel="µg/m³",
    title="(d) Nitrate")
plot(p6a, p6b, p6c, p6d, layout=(2,2), size=(800,600),
    plot_title="Figure 6: Urban (Case 3, Sulfate Rich)")
```

### Figure 7: Marine Aerosol — Sulfate Poor, Crustal Rich (Case 12, Table 8)

Aerosol composition vs. RH for a marine case (R₁ > 2, R₂ > 2):
Na=3.0, H₂SO₄=3.0, NH₃=0.020, HNO₃=2.0, HCl=3.121, Ca=0.36, K=0.45, Mg=0.13 µg/m³.
Panels show (a) aerosol water, (b) K⁺, (c) Na⁺ (metastable — no NaCl solid), and
(d) Mg²⁺. The original Fig. 7(c) shows NaCl(s), which is zero in the metastable state.

```@example iso2
# Case 12 (Marine) from Table 8
RH7, W7, c7, g7 = solve_iso2_sweep(compiled,
    3.0e-6/22.99, 3.0e-6/98.08, 0.020e-6/17.03, 2.0e-6/63.01, 3.121e-6/36.46,
    0.36e-6/40.08, 0.45e-6/39.10, 0.13e-6/24.31)

p7a = plot(RH7, W7 .* 1e9, label="Metastable", xlabel="RH", ylabel="µg/m³",
    title="(a) Aerosol Water", legend=:topleft)
p7b = plot(RH7, c7[:K] .* 39.10e6, label="K⁺", xlabel="RH", ylabel="µg/m³",
    title="(b) Potassium")
p7c = plot(RH7, c7[:Na] .* 22.99e6, label="Na⁺ (aq)", xlabel="RH", ylabel="µg/m³",
    title="(c) Na⁺ (no solids in metastable)")
p7d = plot(RH7, c7[:Mg] .* 24.31e6, label="Mg²⁺", xlabel="RH", ylabel="µg/m³",
    title="(d) Magnesium")
plot(p7a, p7b, p7c, p7d, layout=(2,2), size=(800,600),
    plot_title="Figure 7: Marine (Case 12, Sulfate Poor)")
```

### Figure 8: Non-urban Continental — Sulfate Poor (Case 5, Table 8)

Aerosol composition vs. RH for a non-urban continental case (R₁ > 2, R₂ < 2):
Na=0.2, H₂SO₄=2.0, NH₃=8.0, HNO₃=12.0, HCl=0.2, Ca=0.12, K=0.18, Mg=0 µg/m³.
Panels show (a) aerosol water, (b) NO₃⁻, and (c) NH₄⁺.

```@example iso2
# Case 5 (Non-urban Continental) from Table 8
RH8, W8, c8, g8 = solve_iso2_sweep(compiled,
    0.2e-6/22.99, 2.0e-6/98.08, 8.0e-6/17.03, 12.0e-6/63.01, 0.2e-6/36.46,
    0.12e-6/40.08, 0.18e-6/39.10, 0.0)

p8a = plot(RH8, W8 .* 1e9, label="Metastable", xlabel="RH", ylabel="µg/m³",
    title="(a) Aerosol Water", legend=:topleft)
p8b = plot(RH8, c8[:NO3] .* 62.00e6, label="NO₃⁻", xlabel="RH", ylabel="µg/m³",
    title="(b) Nitrate")
p8c = plot(RH8, c8[:NH4] .* 18.04e6, label="NH₄⁺", xlabel="RH", ylabel="µg/m³",
    title="(c) Ammonium")
plot(p8a, p8b, p8c, layout=(1,3), size=(900,350),
    plot_title="Figure 8: Non-urban Continental (Case 5, Sulfate Poor)")
```

### Figure 9: Remote Continental — Near Neutral (Case 13, Table 8)

Aerosol composition vs. RH for a remote continental case (R₁ ≈ 2):
Na=0, H₂SO₄=10.0, NH₃=4.25, HNO₃=0.145, HCl=0, Ca=0.08, K=0.09, Mg=0 µg/m³.
Panels show (a) aerosol water and (b) K⁺. The metastable branch is shown; set `stable=1`
to compute the stable (solid-forming) branch.

```@example iso2
# Case 13 (Remote Continental) from Table 8
RH9, W9, c9, g9 = solve_iso2_sweep(compiled,
    0.0, 10.0e-6/98.08, 4.25e-6/17.03, 0.145e-6/63.01, 0.0,
    0.08e-6/40.08, 0.09e-6/39.10, 0.0)

p9a = plot(RH9, W9 .* 1e9, label="Metastable", xlabel="RH", ylabel="µg/m³",
    title="(a) Aerosol Water", legend=:topleft)
p9b = plot(RH9, c9[:K] .* 39.10e6, label="K⁺", xlabel="RH", ylabel="µg/m³",
    title="(b) Potassium")
plot(p9a, p9b, layout=(1,2), size=(700,350),
    plot_title="Figure 9: Remote Continental (Case 13, R₁ ≈ 2)")
```
