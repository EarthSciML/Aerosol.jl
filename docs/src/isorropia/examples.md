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


## Helper functions

Then, we need to create some helper functions:

Function to run a sweep through a range of relative humidities:
```@example 1
function run_rh_sweep(sys, RHs, ics; mstable=0)
    defaults = ModelingToolkit.get_defaults(sys)
    u₀ = Dict{Any,Float64}([s => 1.0e-15 for s ∈ keys(model.mw)])
    for k ∈ keys(ics)
        u₀[k] = ics[k] / 1e6 / model.mw[k] # ug/m3 / (1e6 ug/g) / g/mol = mol/m3
    end
    u₀[:H_aq] = 2 * u₀[:SO4_aq]
    p = Dict{Any,Float64}([p => defaults[p] for p ∈ parameters(sys)])
    p[sys.deliquescence.metastable] = mstable

    sols = []
    for rh in RHs
        p[sys.RH] = rh
        local prob = ODEProblem(sys, u₀, (0, 5.0), p)
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
    p1 = plot(RHs, [sols[i][sys.W][end] * 1e9 for i ∈ 1:length(RHs)], #ylim=(0, 50),
    ylabel="H2O (ug/m3)", xlabel="Relative humidity (%)", label=:none)
    ps = []
    for v in plotvars
        push!(ps, plot(RHs, [sols[i][v][end] * 1e6 * model.mw[Symbolics.tosymbol(v, escape=false)] for i ∈ 1:length(RHs)],
            ylabel="$v (ug/m3)", xlabel="Relative humidity (%)", label=:none))
    end
    plot(p1, ps...)
end
nothing # hide
```

Functions to make mass balance plots:
```@example 1
varfromsym(sym) = states(sys)[findfirst((v) -> isequal(Symbolics.tosymbol(v, escape=false), sym), states(sys))]
function plot_mass(RHs, molecs, u₀, sols, title; kwargs...)
    y₀ = zeros(length(RHs))
    y = zeros(length(RHs), length(molecs))
    lab = []
    for (j, (s, x)) in enumerate(molecs)
        y₀ .+= s * u₀[x]
        y[:, j] = s .* [sols[i][varfromsym(x)][end] for i ∈ 1:length(RHs)]
        push!(lab, string(x))
    end
    p1 = plot(ylabel="$(title) (mol/m3)", xlabel="Relative humidity (%)"; kwargs...)
    areaplot!(p1, RHs, y, label=permutedims(lab))
    plot!(p1, RHs, y₀, label="u₀ $(title)", color=:black, linewidth=2)
end

function plot_all_masses(RHs, u₀, sols)
    plot(
        [plot_mass(RHs, model.molecs[molec], u₀, sols, molec) for molec in keys(model.molecs)]...,
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
ics = Dict([:Na_aq => 0, :SO4_aq => 10, :NH3_g => 3.4, :HNO3_g => 2, :HCl_g => 0,
    :Ca_aq => 0.4, :K_aq => 0.33, :Mg_aq => 1e-20]) # ug/m3
u₀, sols = run_rh_sweep(sys, RHs, ics);
nothing # hide
```

```@example 1
plot_rh_sweep(RHs, u₀, sols, [sys.K_aq, sys.NH4_aq, sys.NO3_aq])
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
ics = Dict([:Na_aq => 3, :SO4_aq => 3, :NH3_g => 0.02, :HNO3_g => 2, :HCl_g => 3.121,
    :Ca_aq => 0.360, :K_aq => 0.450, :Mg_aq => 0.130]) # ug/m3
u₀, sols = run_rh_sweep(sys, RHs, ics);
nothing # hide
```

```@example 1
plot_rh_sweep(RHs, u₀, sols, [sys.K_aq, sys.NH4_aq, sys.NO3_aq])
```

```@example 1
plot_all_masses(RHs, u₀, sols)
```

## Non-urban continental aerosol

Here's our representation of Figure 8 from the paper, representing a non-urban continental aerosol case:

```@example 1
ics = Dict([:Na_aq => 0.2, :SO4_aq => 2.0, :NH3_g => 8.0, :HNO3_g => 12, :HCl_g => 0.2,
    :Ca_aq => 0.120, :K_aq => 0.180, :Mg_aq => 0.000]) # ug/m3
u₀, sols = run_rh_sweep(sys, RHs, ics);
nothing # hide
```

```@example 1
plot_rh_sweep(RHs, u₀, sols, [sys.K_aq, sys.NH4_aq, sys.NO3_aq])
```

```@example 1
plot_all_masses(RHs, u₀, sols)
```

## Stable and metastable solutions

Finally, here is our representation of Figure 9 from the paper, showing metastable behavior (where RH is decreasing and salts become supersaturated) vs. stable behavior (where solids can form):

```@example 1
ics = Dict([:Na_aq => 0.0, :SO4_aq => 10.0, :NH3_g => 4.250, :HNO3_g => 0.145, :HCl_g => 0.0,
    :Ca_aq => 0.080, :K_aq => 0.090, :Mg_aq => 0.000]) # ug/m3
u₀1, sols1 = run_rh_sweep(sys, RHs, ics, mstable=0);
u₀2, sols2 = run_rh_sweep(sys, RHs, ics, mstable=1);
nothing # hide
```

```@example 1
p1 = plot(RHs, [sols1[i][sys.W][end] * 1e9 for i ∈ 1:length(RHs)], ylim=(0, 50),
    ylabel="H2O (ug/m3)", xlabel="Relative humidity (%)", label="Stable")
plot!(p1, RHs, [sols2[i][sys.W][end] * 1e9 for i ∈ 1:length(RHs)], ylim=(0, 50), label="Metastable")
p = plot(RHs, [sols1[i][sys.K_aq][end] * 1e6 * model.mw[:K_aq] for i ∈ 1:length(RHs)],
    ylabel="K_aq (ug/m3)", xlabel="Relative humidity (%)", label="Stable")
plot!(p, RHs, [sols2[i][sys.K_aq][end] * 1e6 * model.mw[:K_aq] for i ∈ 1:length(RHs)], label="Metastable")
plot(p1, p, size=(600, 300))
```
