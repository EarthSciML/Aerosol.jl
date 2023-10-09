# [Isorropia: Thermodynamic equilibrium model for K, Ca2, Mg, NH4, Na, SO4, NO3, Cl, and H2O aerosols](@id isorropia)

## Model

This is an implementation of the ISORROPIA II model for thermodynamic equilibrium between gases, inorganic aerosols, solids, and water, as described in Fountoukis and Nenes (2007):

> [Fountoukis, C. and Nenes, A., 2007. ISORROPIA II: a computationally efficient thermodynamic equilibrium model for K, Ca, Mg, NH4, Na, SO4, NO3, Cl, and H2O aerosols. *Atmospheric Chemistry and Physics*, 7(17), pp.4639-4659.](https://doi.org/10.5194/acp-7-4639-2007)

We can create an instance of the model in the following manner:

```@example 1
using Aerosol
using EarthSciMLBase
using ModelingToolkit, Catalyst
using Unitful
using DataFrames

@variables t [unit = u"s", description = "Time"]

model = Isorropia(t);
nothing # hide
```

To get a sense of the complexity involved, we can view a graph of the reaction network:

```@example 1
Graph(model.rxn_sys)
```

Before we run any simulations with the model we need to convert it into a system of differential equations.
Below, we visualize just the first three of them:
```@example 1
sys = structural_simplify(get_mtk(model))
equations(sys)[1:3]
```
## State variables
This system of equations contains the following state variables, which are the variables that will be solved for:

```@example 1
function vars2dataframe(vars; include_defaults=true)
    df = DataFrame(
        :Name => [string(Symbolics.tosymbol(v, escape=false)) for v ∈ vars],
        :Units => [ModelingToolkit.get_unit(v) for v ∈ vars],
        :Description => [ModelingToolkit.getdescription(v) for v ∈ vars])
    if include_defaults
        df.Default = [ModelingToolkit.getdefault(v) for v ∈ vars]
    end
    df
end

vars2dataframe(states(sys))
```
## Parameters

The model also has the following parameters, which are variables that are not solved for but can vary 
from simulation to simulation (there are are also constants, which we are filtering out):

```@example 1
vars2dataframe(parameters(sys)[[!ModelingToolkit.isconstant(p) for p in parameters(sys)]])
```
## Observed variables

Finally, the model has the following observed variables, which are variables that can be solved for by the system
but are not strictly necessary to specify the system state:

```@example 1
vars2dataframe([eq.lhs for eq in observed(sys)]; include_defaults=false)
```
