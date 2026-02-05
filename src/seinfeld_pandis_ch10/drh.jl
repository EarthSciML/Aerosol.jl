"""
Deliquescence Relative Humidity (DRH) equations from Seinfeld & Pandis (2006) Chapter 10.

DRH is the relative humidity at which a crystalline salt spontaneously absorbs water
and forms a saturated aqueous solution.
"""

# DRH values at 298 K from Table 10.1 (values in decimal form)
const DRH_298 = Dict{Symbol, Float64}(
    :KCl => 0.842,
    :Na2SO4 => 0.842,
    :NH4Cl => 0.800,
    :NH42SO4 => 0.799,  # (NH4)2SO4
    :NaCl => 0.753,
    :NaNO3 => 0.743,
    :NH43HSO42 => 0.690,  # (NH4)3H(SO4)2
    :NH4NO3 => 0.618,
    :NaHSO4 => 0.520,
    :NH4HSO4 => 0.400
)

# Enthalpy of solution at 298 K from Table 10.3 (kJ/mol)
const DELTA_HS_298 = Dict{Symbol, Float64}(
    :NH42SO4 => 6.32,    # (NH4)2SO4
    :Na2SO4 => -9.76,
    :NaNO3 => 13.24,
    :NH4NO3 => 16.27,
    :KCl => 15.34,
    :NaCl => 1.88
)

# Solubility parameter n(T) = A + B*T + C*T^2 from Table 10.2
# Returns moles salt per mole water at saturation
const SOLUBILITY_PARAMS = Dict{Symbol, NTuple{3, Float64}}(
    :NH42SO4 => (0.1149, -4.489e-4, 1.385e-6),
    :Na2SO4 => (0.3754, -1.763e-3, 2.424e-6),
    :NaNO3 => (0.1868, -1.677e-3, 5.714e-6),
    :NH4NO3 => (4.298, -3.623e-2, 7.853e-5),
    :KCl => (-0.2368, 1.453e-3, -1.238e-6),
    :NaCl => (0.1805, -5.310e-4, 9.965e-7)
)

"""
    DRHTemperature(; name=:DRHTemperature, salt=:NH4NO3)

A component that calculates the temperature-dependent deliquescence relative humidity.

Implements Equation 10.72 from Seinfeld & Pandis (2006):

```math
DRH(T) = DRH(298) \\exp\\left\\{\\frac{\\Delta H_s}{R}\\left[A\\left(\\frac{1}{T} - \\frac{1}{298}\\right) - B\\ln\\frac{T}{298} - C(T-298)\\right]\\right\\}
```

where A, B, C are the solubility polynomial coefficients from Table 10.2.

# Parameters

  - `T`: Temperature (K)
  - `DRH_298K`: DRH at 298 K (dimensionless, 0-1)
  - `ΔH_s`: Enthalpy of solution (J/mol)
  - `A_sol`: Solubility coefficient A (dimensionless)
  - `B_sol`: Solubility coefficient B (K⁻¹)
  - `C_sol`: Solubility coefficient C (K⁻²)

# Variables

  - `DRH`: Deliquescence relative humidity at temperature T (dimensionless, 0-1)

# Example

```julia
using ModelingToolkit, Aerosol
@mtkbuild sys = DRHTemperature(salt = :NH4NO3)
```
"""
@component function DRHTemperature(; name = :DRHTemperature, salt::Symbol = :NH4NO3)
    @constants begin
        R_gas = 8.314, [description = "Universal gas constant", unit = u"J/mol/K"]
        T_ref = 298.0, [description = "Reference temperature", unit = u"K"]
    end

    # Get salt-specific values from tables
    drh_ref = get(DRH_298, salt, 0.75)
    delta_hs = get(DELTA_HS_298, salt, 0.0) * 1000.0  # Convert kJ to J

    # Get solubility polynomial coefficients from Table 10.2
    sol_params = get(SOLUBILITY_PARAMS, salt, (0.1, 0.0, 0.0))

    @parameters begin
        T = 298.0, [description = "Temperature", unit = u"K"]
        DRH_298K = drh_ref, [description = "DRH at 298 K (dimensionless)", unit = u"1"]
        ΔH_s = delta_hs, [description = "Enthalpy of solution", unit = u"J/mol"]
        A_sol = sol_params[1],
        [description = "Solubility coefficient A (dimensionless)", unit = u"1"]
        B_sol = sol_params[2], [description = "Solubility coefficient B", unit = u"K^-1"]
        C_sol = sol_params[3], [description = "Solubility coefficient C", unit = u"K^-2"]
    end

    @variables begin
        DRH(t),
        [description = "Deliquescence RH at temperature T (dimensionless)", unit = u"1"]
    end

    eqs = [
    # Eq. 10.72: Full form with A, B, C solubility coefficients
        DRH ~
        DRH_298K * exp(ΔH_s / R_gas *
            (A_sol * (1/T - 1/T_ref) - B_sol * log(T/T_ref) - C_sol * (T - T_ref))),  # Eq. 10.72
    ]

    return System(eqs, t; name)
end

"""
    drh_temperature(T, salt::Symbol)

Calculate the DRH of a salt at temperature T using the full form of Eq. 10.72.

```math
DRH(T) = DRH(298) \\exp\\left\\{\\frac{\\Delta H_s}{R}\\left[A\\left(\\frac{1}{T} - \\frac{1}{298}\\right) - B\\ln\\frac{T}{298} - C(T-298)\\right]\\right\\}
```

where A, B, C are the solubility polynomial coefficients from Table 10.2.

# Arguments

  - `T`: Temperature (K)
  - `salt`: Salt symbol (:NH4NO3, :NH42SO4, :NaCl, etc.)

# Returns

DRH as a decimal fraction (0-1)

# Example

```julia
drh = drh_temperature(280.0, :NH4NO3)  # Returns DRH of NH4NO3 at 280 K
```
"""
function drh_temperature(T::Real, salt::Symbol)
    R_gas = 8.314  # J/(mol*K)
    T_ref = 298.0

    drh_ref = get(DRH_298, salt, 0.75)
    delta_hs = get(DELTA_HS_298, salt, 0.0) * 1000.0  # kJ to J

    A, B, C = get(SOLUBILITY_PARAMS, salt, (0.1, 0.0, 0.0))

    # Eq. 10.72: Full form with A, B, C solubility coefficients
    return drh_ref * exp(delta_hs / R_gas *
               (A * (1/T - 1/T_ref) - B * log(T/T_ref) - C * (T - T_ref)))
end

"""
    nh4no3_drh(T)

Calculate the DRH of NH4NO3 using the empirical correlation (Eq. 10.88).

```math
\\ln(DRH) = \\frac{723.7}{T} + 1.6954
```

# Arguments

  - `T`: Temperature (K)

# Returns

DRH as a decimal fraction (0-1)
"""
function nh4no3_drh(T::Real)
    return exp(723.7 / T + 1.6954) / 100.0
end
