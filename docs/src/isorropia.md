# ISORROPIA II

## Overview

ISORROPIA II is a computationally efficient thermodynamic equilibrium model for K⁺–Ca²⁺–Mg²⁺–NH₄⁺–Na⁺–SO₄²⁻–NO₃⁻–Cl⁻–H₂O aerosol systems. The model calculates the equilibrium partitioning between gas, liquid, and solid phases for multicomponent inorganic aerosols under given temperature, relative humidity, and total species concentrations.

**Reference**: Fountoukis, C. and Nenes, A., 2007. ISORROPIA II: a computationally efficient thermodynamic equilibrium model for K⁺–Ca²⁺–Mg²⁺–NH₄⁺–Na⁺–SO₄²⁻–NO₃⁻–Cl⁻–H₂O aerosols. *Atmospheric Chemistry and Physics*, 7(17), pp.4639-4659. [https://doi.org/10.5194/acp-7-4639-2007](https://doi.org/10.5194/acp-7-4639-2007)

See the [`Isorropia`](@ref) API documentation for constructor details.

## Implementation

The ISORROPIA II model is implemented using ModelingToolkit.jl and includes:

- **Thermodynamic equilibrium constants** for 27 chemical reactions (Table 2 of reference paper)
- **Activity coefficient calculations** using multicomponent Debye-Hückel and Kusik-Meissner formulations
- **Aerosol type classification** based on sulfate ratio and crustal species content (Section 3.1)
- **Phase partitioning** between gas, aqueous, and solid phases
- **Temperature dependence** for all equilibrium constants (Equation 5)

### State Variables

```@example isorropia_analysis
using Aerosol, ModelingToolkit, DataFrames, Plots
using OrdinaryDiffEq, NonlinearSolve
using DynamicQuantities, Symbolics

@named sys = Isorropia()

vars = unknowns(sys)
println("Number of state variables: ", length(vars))
# Show a summary of key variable categories
var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in vars]
categories = Dict(
    "Species totals" => count(n -> contains(n, ".total"), var_names),
    "Gas phase" => count(n -> startswith(n, "g₊"), var_names),
    "Aqueous salts" => count(n -> contains(n, ".logm") || contains(n, ".loga"), var_names),
    "Type classification" => count(n -> startswith(n, "type") || startswith(n, "R_"), var_names),
    "Other" => length(vars) - count(n -> contains(n, ".total") || startswith(n, "g₊") || contains(n, ".logm") || contains(n, ".loga") || startswith(n, "type") || startswith(n, "R_"), var_names),
)
DataFrame(:Category => collect(keys(categories)), :Count => collect(values(categories)))
```

## Analysis

### Equilibrium Constants (Table 2)

The equilibrium constants at the reference temperature T₀ = 298.15 K are given in Table 2 of the paper. The following table verifies that the implementation matches the K⁰ values from Table 2.

```@example isorropia_analysis
@named eq = ISORROPIA.EquilibriumConstants()
compiled_eq = mtkcompile(eq)

# K⁰ values from Table 2 of Fountoukis and Nenes (2007)
table2_K0 = [
    (:r1,  "Ca(NO₃)₂(s) ⇌ Ca²⁺ + 2NO₃⁻",      6.067e5),
    (:r2,  "CaCl₂(s) ⇌ Ca²⁺ + 2Cl⁻",            7.974e11),
    (:r3,  "CaSO₄(s) ⇌ Ca²⁺ + SO₄²⁻",           4.319e-5),
    (:r4,  "K₂SO₄(s) ⇌ 2K⁺ + SO₄²⁻",            1.569e-2),
    (:r5,  "KHSO₄(s) ⇌ K⁺ + HSO₄⁻",             24.016),
    (:r6,  "KNO₃(s) ⇌ K⁺ + NO₃⁻",               0.872),
    (:r7,  "KCl(s) ⇌ K⁺ + Cl⁻",                  8.68),
    (:r8,  "MgSO₄(s) ⇌ Mg²⁺ + SO₄²⁻",           1.079e5),
    (:r9,  "Mg(NO₃)₂(s) ⇌ Mg²⁺ + 2NO₃⁻",        2.507e15),
    (:r10, "MgCl₂(s) ⇌ Mg²⁺ + 2Cl⁻",            9.557e21),
    (:r11, "HSO₄⁻ ⇌ H⁺ + SO₄²⁻",                1.015e-2),
    (:r12, "NH₃(g) ⇌ NH₃(aq)",                   5.764e1),
    (:r13, "NH₃(aq) + H₂O ⇌ NH₄⁺ + OH⁻",        1.805e-5),
    (:r14, "HNO₃(g) ⇌ H⁺ + NO₃⁻",               2.511e6),
    (:r16, "HCl(g) ⇌ H⁺ + Cl⁻",                  1.971e6),
    (:r18, "H₂O ⇌ H⁺ + OH⁻",                     1.01e-14),
    (:r19, "Na₂SO₄(s) ⇌ 2Na⁺ + SO₄²⁻",          4.799e-1),
    (:r20, "(NH₄)₂SO₄(s) ⇌ 2NH₄⁺ + SO₄²⁻",      1.817e0),
    (:r22, "NaNO₃(s) ⇌ Na⁺ + NO₃⁻",              1.197e1),
    (:r26, "NH₄HSO₄(s) ⇌ NH₄⁺ + HSO₄⁻",         1.382e2),
]

prob = NonlinearProblem(compiled_eq, [], [compiled_eq.T => 298.15])
sol = solve(prob)

results = DataFrame(
    :Reaction => [r[2] for r in table2_K0],
    :Paper_K0 => [r[3] for r in table2_K0],
    :Computed_K0 => [exp(sol[getproperty(compiled_eq, r[1]).logK_eq]) for r in table2_K0],
)
results.Relative_Error = abs.(results.Computed_K0 .- results.Paper_K0) ./ results.Paper_K0
results
```

### Temperature Dependence of Equilibrium Constants (Equation 5)

The temperature dependence of each equilibrium constant follows Equation 5 of the paper:

``\ln K(T) = \ln K^0 + H^+ \left(\frac{T_0}{T} - 1\right) + C^+ \left(1 + \ln\frac{T_0}{T} - \frac{T_0}{T}\right)``

The following plot shows how selected equilibrium constants vary with temperature over the atmospherically relevant range.

```@example isorropia_analysis
T_range = 250.0:2.0:320.0
reactions_to_plot = [
    (:r11, "HSO₄⁻ dissociation"),
    (:r12, "NH₃(g) ⇌ NH₃(aq)"),
    (:r14, "HNO₃(g) ⇌ H⁺ + NO₃⁻"),
    (:r16, "HCl(g) ⇌ H⁺ + Cl⁻"),
    (:r20, "(NH₄)₂SO₄ dissolution"),
]

p = plot(xlabel = "Temperature (K)", ylabel = "ln K",
    title = "Temperature Dependence of Equilibrium Constants (Eq. 5)",
    legend = :outertopright, size = (800, 400))
for (rname, label) in reactions_to_plot
    logK_vals = Float64[]
    for T in T_range
        prob_T = NonlinearProblem(compiled_eq, [], [compiled_eq.T => T])
        sol_T = solve(prob_T)
        push!(logK_vals, sol_T[getproperty(compiled_eq, rname).logK_eq])
    end
    plot!(p, T_range, logK_vals, label = label, linewidth = 2)
end
p
```

### Deliquescence Relative Humidity vs Temperature (Figure 2)

Figure 2 of the paper shows the deliquescence relative humidity (DRH) as a function of temperature for all ISORROPIA salts. The DRH determines the relative humidity above which a solid salt fully dissolves into the aqueous phase. The temperature dependence is given by:

``\text{DRH}(T) = \text{DRH}(T_0) \exp\left[-l_t \left(\frac{1}{T} - \frac{1}{T_0}\right)\right]``

where ``l_t = -\frac{18}{1000 R} L_s m_s`` (see text after Equation 22 in the paper).

```@example isorropia_analysis
# DRH data from Table 4 of Fountoukis and Nenes (2007)
# (salt name, DRH at 298.15K, l_t parameter)
salts_drh = [
    ("NaCl", 0.7528, 25.0),
    ("Na₂SO₄", 0.93, 80.0),
    ("NaNO₃", 0.7379, 304.0),
    ("NaHSO₄", 0.52, -45.0),
    ("(NH₄)₂SO₄", 0.7997, 80.0),
    ("NH₄NO₃", 0.6183, 852.0),
    ("NH₄Cl", 0.771, 239.0),
    ("NH₄HSO₄", 0.40, 384.0),
    ("(NH₄)₃H(SO₄)₂", 0.69, 186.0),
    ("K₂SO₄", 0.9751, 35.6),
    ("KCl", 0.8426, 158.9),
    ("MgCl₂", 0.3284, 42.23),
    ("Mg(NO₃)₂", 0.54, 230.2),
    ("MgSO₄", 0.8613, -714.5),
    ("Ca(NO₃)₂", 0.4906, 509.4),
    ("CaCl₂", 0.283, 551.1),
]

T_range_drh = 250.0:1.0:320.0
T0 = 298.15

p = plot(xlabel = "Temperature (K)", ylabel = "DRH",
    title = "DRH vs Temperature (cf. Figure 2)",
    legend = :outertopright, size = (900, 500), ylim = (0, 1))

for (name, drh0, lt) in salts_drh
    if isnan(lt)
        # DRH doesn't vary with temperature
        plot!(p, T_range_drh, fill(drh0, length(T_range_drh)), label = name, linewidth = 1.5)
    else
        drh_vals = [drh0 * exp(-lt * (1/T - 1/T0)) for T in T_range_drh]
        plot!(p, T_range_drh, drh_vals, label = name, linewidth = 1.5)
    end
end
p
```

### Model Performance Validation

The following analysis runs the full ISORROPIA model and examines its behavior under varying relative humidity conditions. This demonstrates the gas-aerosol partitioning that is the core function of the model.

#### RH Sensitivity: Aerosol Water Content

As relative humidity increases, aerosol water content increases because more salts deliquesce and take up water. This is a fundamental behavior described in Section 2.3 of the paper.

```@example isorropia_analysis
compiled_sys = mtkcompile(sys)

RH_values = 0.1:0.05:0.95
water_content = Float64[]
gas_NH3 = Float64[]
gas_HNO3 = Float64[]
gas_HCl = Float64[]

for rh in RH_values
    prob = ODEProblem(
        compiled_sys, [], (0.0, 10.0),
        [compiled_sys.RH => rh],
    )
    sol = solve(prob, Rosenbrock23())
    push!(water_content, sol[compiled_sys.aq.W][end])
    push!(gas_NH3, sol[compiled_sys.g.NH3.M][end])
    push!(gas_HNO3, sol[compiled_sys.g.HNO3.M][end])
    push!(gas_HCl, sol[compiled_sys.g.HCl.M][end])
end

p1 = plot(collect(RH_values), water_content .* 1e6,
    xlabel = "Relative Humidity", ylabel = "Aerosol Water (μg/m³ × ρ_w)",
    title = "Aerosol Water Content vs RH",
    linewidth = 2, label = "W", yscale = :log10)

p1
```

#### RH Sensitivity: Gas-Phase Concentrations

As RH increases, volatile species (NH₃, HNO₃, HCl) partition more strongly into the aerosol aqueous phase, reducing their gas-phase concentrations.

```@example isorropia_analysis
p2 = plot(xlabel = "Relative Humidity", ylabel = "Gas Concentration (mol/m³)",
    title = "Gas-Phase Concentrations vs RH",
    legend = :topright, size = (700, 400))
plot!(p2, collect(RH_values), gas_NH3, label = "NH₃(g)", linewidth = 2)
plot!(p2, collect(RH_values), gas_HNO3, label = "HNO₃(g)", linewidth = 2)
plot!(p2, collect(RH_values), gas_HCl, label = "HCl(g)", linewidth = 2)
p2
```

#### Mass Conservation

A key property of the model is that total species concentrations are conserved during the equilibrium calculation (enforced by D(total) ~ 0 for all species). The following verifies this.

```@example isorropia_analysis
prob = ODEProblem(
    compiled_sys, [], (0.0, 10.0),
    [compiled_sys.RH => 0.7],
)
sol = solve(prob, Rosenbrock23())

species = [
    ("NH₃+NH₄⁺", compiled_sys.NH.total),
    ("Na⁺", compiled_sys.Na.total),
    ("Ca²⁺", compiled_sys.Ca.total),
    ("K⁺", compiled_sys.K.total),
    ("Mg²⁺", compiled_sys.Mg.total),
    ("Cl⁻", compiled_sys.Cl.total),
    ("NO₃⁻", compiled_sys.NO3.total),
    ("SO₄²⁻", compiled_sys.SO4.total),
]

DataFrame(
    :Species => [s[1] for s in species],
    :Initial => [sol[s[2]][1] for s in species],
    :Final => [sol[s[2]][end] for s in species],
    :Relative_Change => [abs(sol[s[2]][end] - sol[s[2]][1]) / max(abs(sol[s[2]][1]), 1e-30) for s in species],
)
```
