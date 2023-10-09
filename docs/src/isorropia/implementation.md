This implementation of ISORROPIA II is based on the description by Fountoukis and Nenes (2007),
and may differ from the FORTRAN model which the paper describes. Below, we describe the notable 
differences in model specification and solution that we are aware of.

## Problem formulation

The ISORROPIA II FORTRAN model is formulated as a nonlinear system of equations, which is solved numerically.
For example, the following equation represents the equilibrium between ``\mathrm{Ca(NO3)2_{(s)}}`` and 
``( \mathrm{Ca^{2+}_{(aq)}} + \mathrm{2NO3^-_{(aq)}} )``:

```math
K_{eq} = \frac{\left[ \mathrm{Ca^{2+}} \right] \left[ \mathrm{NO3^-} \right]^2 \gamma_{\mathrm{Ca^{2+}}} \gamma_{NO3^-} }{1}
```


where the numerator (``\left[ \mathrm{Ca^{2+}} \right] \left[ \mathrm{NO3^-} \right]^2 \gamma_{\mathrm{Ca^{2+}}} \gamma_{NO3^-}``) 
represents the activity of the aqueous salt, the denominator (``1``) represents the activity of the solid, 
and ``K_{eq}`` is the equilibrium ratio of the two.

In this implementation, we instead formulate the system of reactions as a mass action problem with forward 
and reverse kinetic reactions that arrive at an equilibrium over time. 
Using this method, the above reaction is represented as:

```math
\mathrm{CaNO32_{s}} \underset{\gamma_p / \gamma_r / K_eq}{\stackrel{1}{\rightleftharpoons}} \mathrm{Ca_{aq}} + 2 \mathrm{NO3_{aq}}
```

where ``( \gamma_p = \gamma_{\mathrm{Ca^{2+}}} \gamma_{NO3^-} )`` and ``( \gamma_r = 1 / \left[ \mathrm{Ca(NO3)2_{(s)}} \right] )``.

The resulting chemical reaction network is converted into a system of ordinary differential equations and then solved.
The advantage of this mass-action formulation is that it is numerically stable even without the extensive 
numerical optimizations that are used in the FORTRAN version of ISORROPIA.
However, we are turning a system that is represented as instantaneously reaching equilibrium into one that
reaches equilibrium over time, so we need to make a choice about the timescale over which we want to reach equilibrium.
We have found that in practice formulating the reaction rates as above results in equilibrium being reached in less than 
10 seconds under typical conditions, which is fast enough compared to many gas-phase reactions that we can consider it 
to be effectively instantaneous.
However, to make the equilibrium occur more quickly or slowly we could simply multiply both forward and reverse reaction 
rate constants by a constant factor.

## Deliquescence

When interpolation between MRDH and DRH relative humidities as in Foungtoukis and Nenes (2007), we take ``\mathrm{RH_{wet}}``
to be equal to the DRH of the species in question, rather than the DRH of the salt with the lowest DRH of the mixture 
under consideration.
(We do this for ease of implementation but would consider changing it if requested.)

## Compositional invariance with RH cycling

This is described in Section 3.2 of Fountoukis and Nenes (2007), but we do not do it 
because it is not immediately clear how to do so in the mass-action framework we use here.