## Running the model

Let's run some simulations!

```@example 1
using Aerosol
using EarthSciMLBase
using ModelingToolkit, Catalyst, DifferentialEquations
using Plots, Latexify, Unitful

@variables t [unit = u"s", description = "Time"]

model = Isorropia(t);

sys = structural_simplify(get_mtk(model))
nothing # hide
```

## Variables

First, we need to extract some variables and parameters from our model to work with:

```@example 1
@unpack Na_aq, SO4_aq, SO4_g, NH3_aq, NH3_g, NO3_aq, Cl_aq, NaCl_s, MgNO32_s, HNO3_g,
    Ca_aq, K_aq, Mg_aq, H_aq, NH4_aq, HCl_g, K2SO4_s, KNO3_s, CaNO32_s, HNO3_g, HNO3_aq,
    KHSO4_s, KCl_s, NH4NO3_s, CaSO4_s, CaCl2_s, MgSO4_s, MgCl2_s, NH4HSO4_s,
    NH42SO4_s, NH43HSO42_s, NH4Cl_s, NaHSO4_s, Na2SO4_s, NaNO3_s, HSO4_aq, HCl_aq,
    RH, metastable, W = sys
nothing # hide
```

## Helper functions

Then, we need to create some helper functions:

Molar masses:
```@example 1 
mw = Dict(Na_aq => 22.989769, SO4_aq => 96.0636, SO4_g => 96.0636, NH3_aq => 17.03052, NH3_g => 17.03052,
    NO3_aq => 62.0049, Cl_aq => 35.453, NaCl_s => 58.44,
    Ca_aq => 40.078, K_aq => 39.0983, Mg_aq => 24.305, H_aq => 1.00784, NH4_aq => 18.04, HCl_g => 36.46,
    K2SO4_s => 174.259, KNO3_s => 101.1032, CaNO32_s => 164.1, HNO3_g => 63.01, HNO3_aq => 63.01,
    KHSO4_s => 136.169, KCl_s => 74.5513, NH4NO3_s => 80.043) # g/mol
nothing # hide
```

Molecules for checking mass balance:
```@example 1
K_molecs = [(1, K_aq), (1, KHSO4_s), (2, K2SO4_s), (1, KNO3_s), (1, KCl_s)]
Ca_molecs = [(1, Ca_aq), (1, CaSO4_s), (1, CaNO32_s), (1, CaCl2_s)]
Mg_molecs = [(1, Mg_aq), (1, MgSO4_s), (1, MgNO32_s), (1, MgCl2_s)]
NH_molecs = [(1, NH4_aq), (1, NH3_aq), (1, NH3_g), (1, NH4HSO4_s),
    (2, NH42SO4_s), (3, NH43HSO42_s), (1, NH4Cl_s), (1, NH4NO3_s)]
Na_molecs = [(1, Na_aq), (1, NaHSO4_s), (2, Na2SO4_s), (1, NaCl_s), (1, NaNO3_s)]
SO4_molecs = [(1, SO4_aq), (1, HSO4_aq), (1, SO4_g),
    (1, KHSO4_s), (1, NaHSO4_s), (1, NH4HSO4_s), (1, CaSO4_s), (1, Na2SO4_s), (1, NH42SO4_s),
    (2, NH43HSO42_s), (1, K2SO4_s), (1, MgSO4_s)]
NO3_molecs = [(1, NO3_aq), (1, HNO3_aq), (1, HNO3_g), (1, NH4NO3_s), (1, NaNO3_s)]
Cl_molecs = [(1, Cl_aq), (1, HCl_aq), (1, HCl_g), (1, NH4Cl_s),
    (1, NaCl_s), (2, CaCl2_s), (1, KCl_s), (2, MgCl2_s)]
H_molecs = [(1, H_aq), (1, HNO3_g), (1, HCl_g), (1, HCl_aq), (1, HSO4_aq)]
nothing # hide
```

Function to run a sweep through a range of relative humidities:
```@example 1
function run_rh_sweep(sys, RHs, ics; mstable=0)
    defaults = ModelingToolkit.get_defaults(sys)
    u₀ = Dict{Any,Float64}([s => 1.0e-15 for s ∈ states(sys)])
    for k ∈ keys(ics)
        u₀[k] = ics[k] / 1e6 / mw[k] # ug/m3 / (1e6 ug/g) / g/mol = mol/m3
    end
    u₀[H_aq] = 2 * u₀[SO4_aq]
    p = Dict{Any,Float64}([p => defaults[p] for p ∈ parameters(sys)])
    p[metastable] = mstable

    sols = []
    for rh in RHs
        p[RH] = rh
        local prob = ODEProblem(sys, u₀, (0, 100.0), p)
        local sol = solve(prob, Rosenbrock23(), abstol=1e-12, reltol=1e-12)
        push!(sols, sol)
    end
    return u₀, sols
end
nothing # hide
```

Function to make plots of our relative humidity sweep:
```@example 1
function plot_rh_sweep(RHs, u₀, sols, plotvars)
    p1 = plot(RHs, [sols[i][W][end] * 1e9 for i ∈ 1:length(RHs)], ylim=(0, 50),
    ylabel="H2O (ug/m3)", xlabel="Relative humidity (%)", label=:none)
    ps = []
    for v in plotvars
        push!(ps, plot(RHs, [sols[i][v][end] * 1e6 * mw[v] for i ∈ 1:length(RHs)],
            ylabel="$v (ug/m3)", xlabel="Relative humidity (%)", label=:none))
    end
    plot(p1, ps...)
end
nothing # hide
```

Functions to make mass balance plots:
```@example 1
function plot_mass(RHs, molecs, u₀, sols, title; kwargs...)
    y₀ = zeros(length(RHs))
    y = zeros(length(RHs), length(molecs))
    lab = []
    for (j, (s, x)) in enumerate(molecs)
        y₀ .+= s * u₀[x]
        y[:, j] = s .* [sols[i][x][end] for i ∈ 1:length(RHs)]
        push!(lab, string(x))
    end
    p1 = plot(ylabel="$(title) (mol/m3)", xlabel="Relative humidity (%)"; kwargs...)
    areaplot!(p1, RHs, y, label=permutedims(lab))
    plot!(p1, RHs, y₀, label="u₀ $(title)", color=:black, linewidth=2)
end

function plot_all_masses(RHs, u₀, sols)
    plot(
        plot_mass(RHs, K_molecs, u₀, sols, "K"),
        plot_mass(RHs, Ca_molecs, u₀, sols, "Ca"),
        plot_mass(RHs, Mg_molecs, u₀, sols, "Mg"),
        plot_mass(RHs, NH_molecs, u₀, sols, "NH"),
        plot_mass(RHs, Na_molecs, u₀, sols, "Na"),
        plot_mass(RHs, SO4_molecs, u₀, sols, "SO4"),
        plot_mass(RHs, NO3_molecs, u₀, sols, "NO3"),
        plot_mass(RHs, Cl_molecs, u₀, sols, "Cl"),
        plot_mass(RHs, H_molecs, u₀, sols, "H"),
        size=(1000, 800)
    )
end
nothing # hide
```

## Urban aerosol

Let's try reproducing Figure 6 by Fountoukis and Nenes (2007).
This case represents an urban aerosol with the following initial conditions:

```@example 1
RHs = [10, 25, 40, 55, 65, 70, 75, 80, 85, 90] ./ 100.0
ics = Dict([Na_aq => 0, SO4_g => 10, NH3_g => 3.4, HNO3_g => 2, HCl_g => 0,
    Ca_aq => 0.4, K_aq => 0.33, Mg_aq => 1e-20]) # ug/m3
u₀, sols = run_rh_sweep(sys, RHs, ics);
nothing # hide
```

```@example 1
plot_rh_sweep(RHs, u₀, sols, [K_aq, NH4_aq, NO3_aq])
```
If you compare it to the paper, it doesn't match exactly, but we're getting there.
Below is a plot of the partitioning between different molecules for each element or ion,
with the black line showing the initial concentration:

```@example 1
plot_all_masses(RHs, u₀, sols)
```

## Marine aerosol

Here's our representation of Figure 7 from the paper, representing a marine aerosol case:

```@example 1
ics = Dict([Na_aq => 3, SO4_g => 3, NH3_g => 0.02, HNO3_g => 2, HCl_g => 3.121,
    Ca_aq => 0.360, K_aq => 0.450, Mg_aq => 0.130]) # ug/m3
u₀, sols = run_rh_sweep(sys, RHs, ics);
nothing # hide
```

```@example 1
plot_rh_sweep(RHs, u₀, sols, [K_aq, NH4_aq, NO3_aq])
```

```@example 1
plot_all_masses(RHs, u₀, sols)
```

## Non-urban continental aerosol

Here's our representation of Figure 8 from the paper, representing a non-urban continental aerosol case:

```@example 1
ics = Dict([Na_aq => 0.2, SO4_g => 2.0, NH3_g => 8.0, HNO3_g => 12, HCl_g => 0.2,
    Ca_aq => 0.120, K_aq => 0.180, Mg_aq => 0.000]) # ug/m3
u₀, sols = run_rh_sweep(sys, RHs, ics);
nothing # hide
```

```@example 1
plot_rh_sweep(RHs, u₀, sols, [K_aq, NH4_aq, NO3_aq])
```

```@example 1
plot_all_masses(RHs, u₀, sols)
```

## Stable and metastable solutions

Finally, here is our representation of Figure 9 from the paper, showing metastable behavior (where RH is decreasing and salts become supersaturated) vs. stable behavior (where solids can form):

```@example 1
ics = Dict([Na_aq => 0.0, SO4_g => 10.0, NH3_g => 4.250, HNO3_g => 0.145, HCl_g => 0.0,
    Ca_aq => 0.080, K_aq => 0.090, Mg_aq => 0.000]) # ug/m3
u₀1, sols1 = run_rh_sweep(sys, RHs, ics, mstable=0);
u₀2, sols2 = run_rh_sweep(sys, RHs, ics, mstable=1);
nothing # hide
```

```@example 1
p1 = plot(RHs, [sols1[i][W][end] * 1e9 for i ∈ 1:length(RHs)], ylim=(0, 50),
    ylabel="H2O (ug/m3)", xlabel="Relative humidity (%)", label="Stable")
plot!(p1, RHs, [sols2[i][W][end] * 1e9 for i ∈ 1:length(RHs)], ylim=(0, 50), label="Metastable")
ps = []
for v in [K_aq]
    p = plot(RHs, [sols1[i][v][end] * 1e6 * mw[v] for i ∈ 1:length(RHs)],
        ylabel="$v (ug/m3)", xlabel="Relative humidity (%)", label="Stable")
    plot!(p, RHs, [sols2[i][v][end] * 1e6 * mw[v] for i ∈ 1:length(RHs)], label="Metastable")
    push!(ps, p)
end
plot(p1, ps..., size=(600, 300))
```
