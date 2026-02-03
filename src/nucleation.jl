export Nucleation, ClassicalNucleationRate, CriticalCluster, WaterProperties

"""
    Nucleation(; name=:Nucleation)

A component representing homogeneous nucleation of a single-component vapor following
classical nucleation theory from Seinfeld & Pandis (2006), Chapter 11.

This component calculates the nucleation rate, critical cluster size, and critical radius
based on the saturation ratio, temperature, and species properties.

# Variables

  - `S`: Saturation ratio (dimensionless), ratio of partial pressure to saturation vapor pressure
  - `J`: Nucleation rate (molecules/(m³·s))
  - `i_star`: Critical cluster size (molecules)
  - `r_star`: Critical cluster radius (m)
  - `ΔG_star`: Free energy barrier (J)
  - `θ`: Dimensionless surface tension parameter

# Parameters

  - `T`: Temperature (K)
  - `p_A`: Partial pressure of condensing species (Pa)
  - `p_A_s`: Saturation vapor pressure (Pa)
  - `σ`: Surface tension (N/m)
  - `v_1`: Molecular volume (m³/molecule)
  - `m_1`: Molecular mass (kg/molecule)

# References

Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 11, pp. 489-536.
"""
@component function Nucleation(; name = :Nucleation)
    @constants begin
        k_B = 1.380649e-23, [description = "Boltzmann constant", unit = u"J/K"]
        π_val = π, [description = "Pi (dimensionless)", unit = u"1"]
    end

    @parameters begin
        T = 293.0, [description = "Temperature", unit = u"K"]
        p_A, [description = "Partial pressure of condensing species", unit = u"Pa"]
        p_A_s, [description = "Saturation vapor pressure", unit = u"Pa"]
        σ, [description = "Surface tension", unit = u"N/m"]
        v_1, [description = "Molecular volume", unit = u"m^3"]
        m_1, [description = "Molecular mass", unit = u"kg"]
    end

    @variables begin
        S(t), [description = "Saturation ratio (dimensionless)", unit = u"1"]
        J(t), [description = "Nucleation rate", unit = u"m^-3*s^-1"]
        i_star(t), [description = "Critical cluster size (dimensionless)", unit = u"1"]
        r_star(t), [description = "Critical cluster radius", unit = u"m"]
        ΔG_star(t), [description = "Free energy barrier", unit = u"J"]
        θ(t), [description = "Dimensionless surface tension (dimensionless)", unit = u"1"]
        N_1(t), [description = "Monomer number concentration", unit = u"m^-3"]
        a_1(t), [description = "Monomer surface area", unit = u"m^2"]
    end

    eqs = [
        # Eq. 11.1: Saturation ratio
        S ~ p_A / p_A_s,

        # Monomer number concentration from ideal gas law: N_1 = p_A / (k_B * T)
        N_1 ~ p_A / (k_B * T),

        # Eq. 11.15: Monomer surface area
        a_1 ~ 4 * π_val * (3 * v_1 / (4 * π_val))^(2 // 3),

        # Eq. 11.27: Dimensionless surface tension parameter
        θ ~ (36 * π_val)^(1 // 3) * v_1^(2 // 3) * σ / (k_B * T),

        # Eq. 11.35: Critical cluster size
        # i* = (2θ / (3 ln S))^3
        i_star ~ (2 * θ / (3 * log(S)))^3,

        # Eq. 11.52: Critical radius
        # r* = 2σv_1 / (kT ln S)
        r_star ~ 2 * σ * v_1 / (k_B * T * log(S)),

        # Eq. 11.53: Free energy barrier
        # ΔG* = (16π/3) * v_1² * σ³ / ((kT ln S)²)
        ΔG_star ~ (16 * π_val / 3) * v_1^2 * σ^3 / ((k_B * T * log(S))^2),

        # Eq. 11.47: Classical nucleation rate
        # J = (2σ/(πm_1))^(1/2) * (v_1 * N_1² / S) * exp(-(16π/3) * v_1² * σ³ / ((kT)³ (ln S)²))
        J ~
        sqrt(2 * σ / (π_val * m_1)) * (v_1 * N_1^2 / S) *
        exp(-(16 * π_val / 3) * v_1^2 * σ^3 / ((k_B * T)^3 * log(S)^2))
    ]

    return System(eqs, t; name)
end

"""
    WaterProperties(; name=:WaterProperties)

A component providing temperature-dependent physical properties of water for nucleation
calculations following the data from Seinfeld & Pandis (2006), Table 11.1 and 11.4.

# Variables

  - `σ`: Surface tension (N/m)
  - `v_1`: Molecular volume (m³/molecule)
  - `m_1`: Molecular mass (kg/molecule)
  - `p_sat`: Saturation vapor pressure (Pa)

# Parameters

  - `T`: Temperature (K)

# References

Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics*, Chapter 11.
"""
@component function WaterProperties(; name = :WaterProperties)
    @constants begin
        # Molecular mass of water
        M_H2O = 18.015e-3, [description = "Molar mass of water", unit = u"kg/mol"]
        N_A = 6.02214076e23, [description = "Avogadro number", unit = u"mol^-1"]
        R = 8.314, [description = "Gas constant", unit = u"J/(mol*K)"]

        # Reference values for surface tension (linear interpolation coefficients)
        # σ ≈ 75.6 - 0.1454*(T-273) mN/m based on Table 11.1 data
        σ_ref = 75.6e-3, [description = "Surface tension at 273 K", unit = u"N/m"]
        σ_slope = -0.1454e-3,
        [description = "Surface tension temperature slope", unit = u"N/(m*K)"]
        T_ref = 273.0, [description = "Reference temperature", unit = u"K"]

        # Molecular volume (approximately constant over this T range)
        # v_1 = 2.99e-29 m³/molecule from Table 11.1
        v_1_const = 2.99e-29, [description = "Molecular volume of water", unit = u"m^3"]

        # Antoine equation constants for water (dimensionless form)
        # ln(p/Pa) = A - B/(C+T/K) where T is in K
        A_antoine = 23.1964,
        [description = "Antoine A coefficient (dimensionless)", unit = u"1"]
        B_antoine = 3816.44,
        [description = "Antoine B coefficient (dimensionless)", unit = u"1"]
        C_antoine = -46.13,
        [description = "Antoine C coefficient (dimensionless)", unit = u"1"]
    end

    @parameters begin
        T = 293.0, [description = "Temperature", unit = u"K"]
        T_unit = 1.0,
        [description = "Unit temperature for conversion (dimensionless)", unit = u"K"]
        Pa_unit = 1.0,
        [description = "Unit pressure for conversion (dimensionless)", unit = u"Pa"]
    end

    @variables begin
        σ(t), [description = "Surface tension", unit = u"N/m"]
        v_1(t), [description = "Molecular volume", unit = u"m^3"]
        m_1(t), [description = "Molecular mass", unit = u"kg"]
        p_sat(t), [description = "Saturation vapor pressure", unit = u"Pa"]
    end

    eqs = [
        # Temperature-dependent surface tension (linear approximation)
        σ ~ σ_ref + σ_slope * (T - T_ref),

        # Molecular volume (approximately constant)
        v_1 ~ v_1_const,

        # Molecular mass
        m_1 ~ M_H2O / N_A,

        # Saturation vapor pressure from Antoine equation
        # Using natural log: ln(p/Pa) = A - B/(C + T/K)
        p_sat ~ Pa_unit * exp(A_antoine - B_antoine / (C_antoine + T / T_unit))
    ]

    return System(eqs, t; name)
end

"""
    CriticalCluster(; name=:CriticalCluster)

A standalone component for computing critical cluster properties without the full
nucleation rate calculation. Useful for parameter sensitivity studies.

This implements equations 11.35, 11.52, and 11.53 from Seinfeld & Pandis (2006).

# Variables

  - `i_star`: Critical cluster size (molecules)
  - `r_star`: Critical cluster radius (m)
  - `ΔG_star`: Free energy barrier (J)
  - `θ`: Dimensionless surface tension parameter

# Parameters

  - `T`: Temperature (K)
  - `S`: Saturation ratio
  - `σ`: Surface tension (N/m)
  - `v_1`: Molecular volume (m³/molecule)
"""
@component function CriticalCluster(; name = :CriticalCluster)
    @constants begin
        k_B = 1.380649e-23, [description = "Boltzmann constant", unit = u"J/K"]
        π_val = π, [description = "Pi (dimensionless)", unit = u"1"]
    end

    @parameters begin
        T = 293.0, [description = "Temperature", unit = u"K"]
        S = 2.0, [description = "Saturation ratio (dimensionless)", unit = u"1"]
        σ, [description = "Surface tension", unit = u"N/m"]
        v_1, [description = "Molecular volume", unit = u"m^3"]
    end

    @variables begin
        i_star(t), [description = "Critical cluster size (dimensionless)", unit = u"1"]
        r_star(t), [description = "Critical cluster radius", unit = u"m"]
        ΔG_star(t), [description = "Free energy barrier", unit = u"J"]
        θ(t), [description = "Dimensionless surface tension (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 11.27: Dimensionless surface tension
        θ ~ (36 * π_val)^(1 // 3) * v_1^(2 // 3) * σ / (k_B * T),

        # Eq. 11.35: Critical cluster size
        i_star ~ (2 * θ / (3 * log(S)))^3,

        # Eq. 11.52: Critical radius
        r_star ~ 2 * σ * v_1 / (k_B * T * log(S)),

        # Eq. 11.53: Free energy barrier
        ΔG_star ~ (16 * π_val / 3) * v_1^2 * σ^3 / ((k_B * T * log(S))^2)
    ]

    return System(eqs, t; name)
end

"""
    ClassicalNucleationRate(; name=:ClassicalNucleationRate)

A simplified component that computes only the classical nucleation rate (Eq. 11.47)
from Seinfeld & Pandis (2006). This is a minimal component for use when critical
cluster properties are not needed.

# Variables

  - `J`: Nucleation rate (molecules/(m³·s))

# Parameters

  - `T`: Temperature (K)
  - `S`: Saturation ratio
  - `N_1`: Monomer number concentration (molecules/m³)
  - `σ`: Surface tension (N/m)
  - `v_1`: Molecular volume (m³/molecule)
  - `m_1`: Molecular mass (kg/molecule)
"""
@component function ClassicalNucleationRate(; name = :ClassicalNucleationRate)
    @constants begin
        k_B = 1.380649e-23, [description = "Boltzmann constant", unit = u"J/K"]
        π_val = π, [description = "Pi (dimensionless)", unit = u"1"]
    end

    @parameters begin
        T = 293.0, [description = "Temperature", unit = u"K"]
        S = 2.0, [description = "Saturation ratio (dimensionless)", unit = u"1"]
        N_1, [description = "Monomer number concentration", unit = u"m^-3"]
        σ, [description = "Surface tension", unit = u"N/m"]
        v_1, [description = "Molecular volume", unit = u"m^3"]
        m_1, [description = "Molecular mass", unit = u"kg"]
    end

    @variables begin
        J(t), [description = "Nucleation rate", unit = u"m^-3*s^-1"]
    end

    eqs = [
    # Eq. 11.47: Classical nucleation rate
        J ~
        sqrt(2 * σ / (π_val * m_1)) * (v_1 * N_1^2 / S) *
        exp(-(16 * π_val / 3) * v_1^2 * σ^3 / ((k_B * T)^3 * log(S)^2)),
    ]

    return System(eqs, t; name)
end
