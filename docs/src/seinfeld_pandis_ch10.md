# Aerosol Thermodynamics (Seinfeld & Pandis Ch. 10)

## Overview

This module implements key aerosol thermodynamics equations from Chapter 10 of Seinfeld & Pandis (2006). These equations govern the phase behavior, hygroscopic growth, and gas-particle partitioning of atmospheric aerosols, which are essential processes for air quality modeling and climate science.

The module provides four ModelingToolkit components:

  - **KelvinEffect**: Vapor pressure enhancement over curved droplet surfaces (Eq. 10.86)
  - **DRHTemperature**: Temperature-dependent deliquescence relative humidity (Eq. 10.72)
  - **ZSRWaterContent**: Aerosol liquid water content via the Zdanovskii-Stokes-Robinson mixing rule (Eq. 10.98)
  - **NH4NO3Equilibrium**: Ammonium nitrate gas-aerosol partitioning (Eqs. 10.88, 10.91)

Additionally, standalone functions are provided for all equilibrium constant calculations from Table 10.7.

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006) "Atmospheric Chemistry and Physics: From Air Pollution to Climate Change", 2nd Edition, John Wiley & Sons, Chapter 10.

```@docs
KelvinEffect
DRHTemperature
ZSRWaterContent
NH4NO3Equilibrium
```

## Implementation

### KelvinEffect Component

The Kelvin effect describes the enhancement of vapor pressure over a curved liquid surface relative to a flat surface (Eq. 10.86):

```math
p_A = p°_A \exp\left(\frac{2\sigma M}{RT \rho_l R_p}\right)
```

#### Variables

```@example seinfeld_pandis_ch10
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Aerosol, Aerosol.SeinfeldPandisCh10

sys = KelvinEffect()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Parameters

```@example seinfeld_pandis_ch10
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

#### Equations

```@example seinfeld_pandis_ch10
equations(sys)
```

### DRHTemperature Component

Calculates the temperature-dependent deliquescence relative humidity using Eq. 10.72:

```math
DRH(T) = DRH(298) \exp\left\{\frac{\Delta H_s}{R}\left[A\left(\frac{1}{T} - \frac{1}{298}\right) - B\ln\frac{T}{298} - C(T-298)\right]\right\}
```

where A, B, C are the solubility polynomial coefficients from Table 10.2.

#### Variables and Parameters

```@example seinfeld_pandis_ch10
sys_drh = DRHTemperature(salt = :NH4NO3)
params_drh = parameters(sys_drh)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_drh],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params_drh],
    :Description => [ModelingToolkit.getdescription(p) for p in params_drh]
)
```

#### Equations

```@example seinfeld_pandis_ch10
equations(sys_drh)
```

### ZSRWaterContent Component

Calculates aerosol liquid water content using the ZSR mixing rule (Eq. 10.98):

```math
W = \sum_i \frac{C_i}{m_{i,0}(a_w)}
```

where ``C_i`` is the concentration of species ``i`` and ``m_{i,0}(a_w)`` is the molality in binary solution at water activity ``a_w = RH`` (Eq. 10.63).

#### Variables and Parameters

```@example seinfeld_pandis_ch10
sys_zsr = ZSRWaterContent(n_species = 2)
vars_zsr = unknowns(sys_zsr)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_zsr],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars_zsr],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_zsr]
)
```

### NH4NO3Equilibrium Component

Implements NH₄NO₃ gas-aerosol equilibrium including:

  - Solid dissociation constant ``K_p`` (Eq. 10.91): ``\ln K_p = 84.6 - 24220/T - 6.1 \ln(T/298)``
  - Empirical DRH (Eq. 10.88): ``\ln(DRH\%) = 723.7/T + 1.6954``
  - Phase determination based on RH vs DRH

#### Equations

```@example seinfeld_pandis_ch10
sys_nh4 = NH4NO3Equilibrium()
equations(sys_nh4)
```

### Equilibrium Constants (Table 10.7)

Temperature-dependent equilibrium constants for 13 reactions follow the general form:

```math
K(T) = K(298) \exp\left\{a\left(\frac{298}{T} - 1\right) + b\left[1 + \ln\frac{298}{T} - \frac{298}{T}\right]\right\}
```

```@docs
get_equilibrium_constant
```

## Analysis

### Kelvin Effect: Saturation Ratio vs. Droplet Radius (cf. Figure 10.12)

The Kelvin effect becomes significant for particles smaller than ~100 nm. This figure shows the saturation ratio as a function of droplet radius for water at 298 K.

```@example seinfeld_pandis_ch10
using Plots

# Saturation ratio vs droplet radius for water (Fig. 10.12)
R_p_range = 10 .^ range(log10(2e-9), log10(1e-6), length = 200)
S_water = [kelvin_saturation_ratio(298.0, R_p, 0.072, 0.018, 1000.0) for R_p in R_p_range]

plot(R_p_range * 1e9, S_water,
    xlabel = "Droplet Radius (nm)", ylabel = "Saturation Ratio S",
    label = "Water (σ=0.072 N/m)",
    xscale = :log10, legend = :topright,
    title = "Kelvin Effect: Vapor Pressure Enhancement",
    linewidth = 2, xlims = (1, 1000), ylims = (1.0, 1.3))

# Add other compounds from Table 10.6 for comparison
# Dioctyl phthalate (DOP)-like organic: lower surface tension
S_organic = [kelvin_saturation_ratio(298.0, R_p, 0.030, 0.104, 900.0) for R_p in R_p_range]
plot!(R_p_range * 1e9, S_organic, label = "Organic (σ=0.030 N/m)",
    linewidth = 2, linestyle = :dash)

# Add reference line at S=1
hline!([1.0], label = "", color = :gray, linestyle = :dot)

savefig("kelvin_effect.svg");
nothing # hide
```

![](kelvin_effect.svg)

### DRH Temperature Dependence (cf. Figure 10.6)

The deliquescence relative humidity varies with temperature. Salts with positive enthalpy of solution (``\Delta H_s > 0``) show decreasing DRH with increasing temperature, while those with negative ``\Delta H_s`` show the opposite trend.

```@example seinfeld_pandis_ch10
T_range = 250.0:1.0:320.0

salts = [:NH4NO3, :NH42SO4, :NaCl, :Na2SO4, :KCl, :NaNO3]
salt_labels = ["NH₄NO₃", "(NH₄)₂SO₄", "NaCl", "Na₂SO₄", "KCl", "NaNO₃"]

p_drh = plot(xlabel = "Temperature (K)", ylabel = "DRH (%)",
    title = "DRH Temperature Dependence (Eq. 10.72)",
    legend = :outertopright, ylims = (40, 100))

for (salt, label) in zip(salts, salt_labels)
    drh_vals = [drh_temperature(T, salt) * 100 for T in T_range]
    plot!(p_drh, T_range, drh_vals, label = label, linewidth = 2)
end

savefig("drh_temperature.svg");
nothing # hide
```

![](drh_temperature.svg)

### NH₄NO₃ Dissociation Constant vs. Temperature (cf. Figure 10.19)

The dissociation constant ``K_p`` of solid NH₄NO₃ increases strongly with temperature, indicating that higher temperatures favor evaporation of NH₄NO₃ to NH₃(g) + HNO₃(g).

```@example seinfeld_pandis_ch10
T_range_kp = 270.0:0.5:320.0
Kp_vals = [nh4no3_Kp(T) for T in T_range_kp]

plot(T_range_kp, Kp_vals,
    xlabel = "Temperature (K)", ylabel = "Kp (ppb²)",
    title = "NH₄NO₃ Dissociation Constant (Eq. 10.91)",
    label = "Kp", linewidth = 2, yscale = :log10,
    legend = :topleft)

savefig("nh4no3_kp.svg");
nothing # hide
```

![](nh4no3_kp.svg)

### NH₄NO₃ Deliquescence RH vs. Temperature (cf. Figure 10.20)

The DRH of NH₄NO₃ decreases with increasing temperature, which determines the phase transition boundary between solid and aqueous NH₄NO₃.

```@example seinfeld_pandis_ch10
T_range_drh = 270.0:0.5:320.0
drh_nh4no3_emp = [nh4no3_drh(T) * 100 for T in T_range_drh]
drh_nh4no3_eq72 = [drh_temperature(T, :NH4NO3) * 100 for T in T_range_drh]

plot(T_range_drh, drh_nh4no3_emp,
    xlabel = "Temperature (K)", ylabel = "DRH (%)",
    title = "NH₄NO₃ DRH vs. Temperature",
    label = "Empirical (Eq. 10.88)", linewidth = 2)
plot!(T_range_drh, drh_nh4no3_eq72,
    label = "General form (Eq. 10.72)", linewidth = 2, linestyle = :dash)

savefig("nh4no3_drh.svg");
nothing # hide
```

![](nh4no3_drh.svg)

### ZSR Water Content vs. Relative Humidity

The ZSR mixing rule predicts aerosol liquid water content as a function of relative humidity for multicomponent aerosols.

```@example seinfeld_pandis_ch10
RH_range = 0.4:0.01:0.95

# Single-component water uptake
W_nh42so4 = [zsr_water_content(RH, Dict(:NH42SO4 => 1e-6)) for RH in RH_range]
W_nh4no3 = [zsr_water_content(RH, Dict(:NH4NO3 => 1e-6)) for RH in RH_range]
W_nacl = [zsr_water_content(RH, Dict(:NaCl => 1e-6)) for RH in RH_range]

# Mixed aerosol
W_mixed = [zsr_water_content(RH, Dict(:NH42SO4 => 0.5e-6, :NH4NO3 => 0.5e-6))
           for RH in RH_range]

plot(RH_range * 100, W_nh42so4 * 1e6,
    xlabel = "Relative Humidity (%)", ylabel = "Water Content (μmol/m³)",
    title = "ZSR Aerosol Water Content (Eq. 10.98)",
    label = "(NH₄)₂SO₄", linewidth = 2)
plot!(RH_range * 100, W_nh4no3 * 1e6, label = "NH₄NO₃", linewidth = 2)
plot!(RH_range * 100, W_nacl * 1e6, label = "NaCl", linewidth = 2)
plot!(RH_range * 100, W_mixed * 1e6, label = "Mixed", linewidth = 2, linestyle = :dash)

savefig("zsr_water.svg");
nothing # hide
```

![](zsr_water.svg)

### Temperature Dependence of Equilibrium Constants (Table 10.7)

Selected equilibrium constants from Table 10.7 showing their temperature dependence.

```@example seinfeld_pandis_ch10
T_range_eq = 260.0:1.0:320.0

reactions = [:HSO4_dissoc, :NaCl_dissoc, :NH42SO4_dissoc, :Na2SO4_dissoc]
labels = ["HSO₄⁻ dissoc.", "NaCl dissoc.", "(NH₄)₂SO₄ dissoc.", "Na₂SO₄ dissoc."]

p_eq = plot(xlabel = "Temperature (K)", ylabel = "K / K(298)",
    title = "Equilibrium Constants (Table 10.7, normalized)",
    legend = :topleft, yscale = :log10)

for (rxn, label) in zip(reactions, labels)
    K_ref = get_equilibrium_constant(rxn, 298.0)
    K_vals = [get_equilibrium_constant(rxn, T) / K_ref for T in T_range_eq]
    plot!(p_eq, T_range_eq, K_vals, label = label, linewidth = 2)
end

hline!([1.0], label = "", color = :gray, linestyle = :dot)

savefig("equilibrium_constants.svg");
nothing # hide
```

![](equilibrium_constants.svg)

### Ionic Strength Fraction (Eq. 10.100)

The ionic strength fraction ``Y`` determines the activity coefficient correction for NH₄NO₃ in the presence of (NH₄)₂SO₄:

```math
Y = \frac{[\text{NH}_4\text{NO}_3]}{[\text{NH}_4\text{NO}_3] + 3[(\text{NH}_4)_2\text{SO}_4]}
```

```@example seinfeld_pandis_ch10
ratios = 0.0:0.01:1.0  # NH4NO3 fraction of total
Y_vals = [ionic_strength_fraction(r, 1-r) for r in ratios]

plot(ratios * 100, Y_vals,
    xlabel = "NH₄NO₃ Mole Fraction (%)", ylabel = "Ionic Strength Fraction Y",
    title = "Ionic Strength Fraction (Eq. 10.100)",
    label = "Y", linewidth = 2, legend = :topleft)

savefig("ionic_strength.svg");
nothing # hide
```

![](ionic_strength.svg)

## Functional API Reference

```@docs
kelvin_saturation_ratio
drh_temperature
nh4no3_drh
zsr_water_content
nh4no3_Kp
nh4no3_K_AN
ionic_strength_fraction
Aerosol.SeinfeldPandisCh10.equilibrium_constant_temperature
Aerosol.SeinfeldPandisCh10.binary_molality_nh42so4
Aerosol.SeinfeldPandisCh10.binary_molality_nh4no3
Aerosol.SeinfeldPandisCh10.binary_molality_nacl
```

## Module Reference

```@docs
Aerosol.SeinfeldPandisCh10
```
