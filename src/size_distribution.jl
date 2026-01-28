export AerosolDistribution
export UrbanAerosol, MarineAerosol, RuralAerosol, RemoteContinentalAerosol
export FreeTroposphereAerosol, PolarAerosol, DesertAerosol

"""
    AerosolDistribution(n_modes=3; name=:AerosolDistribution)

Multi-modal lognormal aerosol size distribution following the parameterization of
Seinfeld and Pandis (2006), Chapter 8, "Properties of the Atmospheric Aerosol",
in *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*,
2nd Edition, John Wiley & Sons.

The model represents the aerosol number, surface area, and volume
distributions as a sum of `n_modes` lognormal modes (Eq. 8.54).
Each mode `i` is characterized by three parameters:
- `N[i]`: total number concentration of mode `i`
- `D_g[i]`: geometric median diameter of mode `i`
- `logσ[i]`: log10 of the geometric standard deviation of mode `i`

The implementation also includes vertical profiles for aerosol mass (Eq. 8.55)
and number concentration (Eq. 8.56).

Default parameter values correspond to the "Urban" distribution from Table 8.3.
"""
@component function AerosolDistribution(n_modes = 3; name = :AerosolDistribution)
    @constants begin
        π_c = π, [description = "Pi", unit = u"1"]
        ln10 = log(10.0), [description = "Natural log of 10", unit = u"1"]
    end

    @parameters begin
        N[1:n_modes], [description = "Number concentration of mode i", unit = u"m^-3"]
        D_g[1:n_modes], [description = "Geometric median diameter of mode i", unit = u"m"]
        logσ[1:n_modes], [description = "Log10 of geometric standard deviation of mode i (dimensionless)", unit = u"1"]
        D_p, [description = "Particle diameter evaluation point", unit = u"m"]
        ρ_p = 1500.0, [description = "Particle density", unit = u"kg/m^3"]
        z = 0.0, [description = "Altitude above surface", unit = u"m"]
        H_p = 1000.0, [description = "Aerosol mass scale height", unit = u"m"]
    end

    @variables begin
        n_N_o(t), [description = "Number distribution dN/d(log D_p)", unit = u"m^-3"]
        n_S_o(t), [description = "Surface area distribution dS/d(log D_p)", unit = u"m^-1"]
        n_V_o(t), [description = "Volume distribution dV/d(log D_p)", unit = u"1"]
        N_t(t), [description = "Total number concentration", unit = u"m^-3"]
        S_t(t), [description = "Total surface area concentration", unit = u"m^-1"]
        V_t(t), [description = "Total volume concentration", unit = u"1"]
        D_s(t), [description = "Surface area median diameter", unit = u"m"]
        D_v(t), [description = "Volume median diameter", unit = u"m"]
        D_bar(t), [description = "Mean diameter", unit = u"m"]
        M_z(t), [description = "Mass concentration vertical profile ratio M(z)/M(0) (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 8.54 — Multi-modal lognormal number distribution dN/d(log D_p)
        n_N_o ~ sum(
            N[i] / (sqrt(2 * π_c) * logσ[i]) *
            exp(-(log10(D_p / D_g[i]))^2 / (2 * logσ[i]^2))
            for i in 1:n_modes
        ),

        # Eq. 8.48/8.50 — Surface area distribution dS/d(log D_p)
        n_S_o ~ sum(
            π_c * N[i] * (D_g[i] * exp(2 * (logσ[i] * ln10)^2))^2 /
            (sqrt(2 * π_c) * logσ[i]) *
            exp(-(log10(D_p / (D_g[i] * exp(2 * (logσ[i] * ln10)^2))))^2 /
                (2 * logσ[i]^2))
            for i in 1:n_modes
        ),

        # Eq. 8.51 — Volume distribution dV/d(log D_p)
        n_V_o ~ sum(
            (π_c / 6) * N[i] * (D_g[i] * exp(3 * (logσ[i] * ln10)^2))^3 /
            (sqrt(2 * π_c) * logσ[i]) *
            exp(-(log10(D_p / (D_g[i] * exp(3 * (logσ[i] * ln10)^2))))^2 /
                (2 * logσ[i]^2))
            for i in 1:n_modes
        ),

        # Total number concentration — sum of all modes
        N_t ~ sum(N[i] for i in 1:n_modes),

        # Total surface area — Eq. 8.5 with Eq. 8.41 (k=2 moment)
        S_t ~ sum(
            π_c * N[i] * D_g[i]^2 * exp(2 * (logσ[i] * ln10)^2)
            for i in 1:n_modes
        ),

        # Total volume — Eq. 8.7 with Eq. 8.41 (k=3 moment)
        V_t ~ sum(
            (π_c / 6) * N[i] * D_g[i]^3 * exp(4.5 * (logσ[i] * ln10)^2)
            for i in 1:n_modes
        ),

        # Eq. 8.49 — Surface area median diameter (first mode)
        D_s ~ D_g[1] * exp(2 * (logσ[1] * ln10)^2),

        # Eq. 8.52 — Volume median diameter (first mode)
        D_v ~ D_g[1] * exp(3 * (logσ[1] * ln10)^2),

        # Eq. 8.39 — Mean diameter (first mode)
        D_bar ~ D_g[1] * exp(0.5 * (logσ[1] * ln10)^2),

        # Eq. 8.55 — Vertical mass profile
        M_z ~ exp(-z / H_p),
    ]

    return System(eqs, t; name)
end

# -------------------------------------------------------------------
# Predefined aerosol distributions from Table 8.3
# (Seinfeld & Pandis, 2006, Chapter 8)
# -------------------------------------------------------------------
# Unit conversions: N from cm^{-3} to m^{-3} (×1e6), D_p from μm to m (×1e-6)

"""
    UrbanAerosol(; name=:UrbanAerosol)

Urban aerosol size distribution from Table 8.3 of Seinfeld and Pandis (2006).

Three lognormal modes with parameters:
- Mode 1: nucleation mode (N=9.93×10⁴ cm⁻³, Dg=0.013 μm, log σ=0.245)
- Mode 2: Aitken mode (N=1.11×10³ cm⁻³, Dg=0.014 μm, log σ=0.666)
- Mode 3: accumulation mode (N=3.64×10⁴ cm⁻³, Dg=0.050 μm, log σ=0.337)
"""
@component function UrbanAerosol(; name = :UrbanAerosol)
    sys = AerosolDistribution(3; name)
    defaults = Dict(
        sys.N[1] => 9.93e4 * 1e6, sys.D_g[1] => 0.013e-6, sys.logσ[1] => 0.245,
        sys.N[2] => 1.11e3 * 1e6, sys.D_g[2] => 0.014e-6, sys.logσ[2] => 0.666,
        sys.N[3] => 3.64e4 * 1e6, sys.D_g[3] => 0.050e-6, sys.logσ[3] => 0.337,
    )
    return System(equations(sys), t; name, defaults,
        systems = ModelingToolkit.get_systems(sys))
end

"""
    MarineAerosol(; name=:MarineAerosol)

Marine aerosol size distribution from Table 8.3 of Seinfeld and Pandis (2006).

Three lognormal modes with parameters:
- Mode 1: (N=133 cm⁻³, Dg=0.008 μm, log σ=0.657)
- Mode 2: (N=66.6 cm⁻³, Dg=0.266 μm, log σ=0.210)
- Mode 3: (N=3.06 cm⁻³, Dg=0.580 μm, log σ=0.396)
"""
@component function MarineAerosol(; name = :MarineAerosol)
    sys = AerosolDistribution(3; name)
    defaults = Dict(
        sys.N[1] => 133.0 * 1e6, sys.D_g[1] => 0.008e-6, sys.logσ[1] => 0.657,
        sys.N[2] => 66.6 * 1e6,  sys.D_g[2] => 0.266e-6, sys.logσ[2] => 0.210,
        sys.N[3] => 3.06 * 1e6,  sys.D_g[3] => 0.580e-6, sys.logσ[3] => 0.396,
    )
    return System(equations(sys), t; name, defaults,
        systems = ModelingToolkit.get_systems(sys))
end

"""
    RuralAerosol(; name=:RuralAerosol)

Rural aerosol size distribution from Table 8.3 of Seinfeld and Pandis (2006).

Three lognormal modes with parameters:
- Mode 1: (N=6.65×10³ cm⁻³, Dg=0.015 μm, log σ=0.225)
- Mode 2: (N=147 cm⁻³, Dg=0.054 μm, log σ=0.557)
- Mode 3: (N=1990 cm⁻³, Dg=0.084 μm, log σ=0.266)
"""
@component function RuralAerosol(; name = :RuralAerosol)
    sys = AerosolDistribution(3; name)
    defaults = Dict(
        sys.N[1] => 6.65e3 * 1e6, sys.D_g[1] => 0.015e-6, sys.logσ[1] => 0.225,
        sys.N[2] => 147.0 * 1e6,  sys.D_g[2] => 0.054e-6, sys.logσ[2] => 0.557,
        sys.N[3] => 1990.0 * 1e6, sys.D_g[3] => 0.084e-6, sys.logσ[3] => 0.266,
    )
    return System(equations(sys), t; name, defaults,
        systems = ModelingToolkit.get_systems(sys))
end

"""
    RemoteContinentalAerosol(; name=:RemoteContinentalAerosol)

Remote continental aerosol size distribution from Table 8.3 of Seinfeld and Pandis (2006).

Three lognormal modes with parameters:
- Mode 1: (N=3200 cm⁻³, Dg=0.020 μm, log σ=0.161)
- Mode 2: (N=2900 cm⁻³, Dg=0.116 μm, log σ=0.217)
- Mode 3: (N=0.300 cm⁻³, Dg=1.800 μm, log σ=0.380)
"""
@component function RemoteContinentalAerosol(; name = :RemoteContinentalAerosol)
    sys = AerosolDistribution(3; name)
    defaults = Dict(
        sys.N[1] => 3200.0 * 1e6, sys.D_g[1] => 0.020e-6, sys.logσ[1] => 0.161,
        sys.N[2] => 2900.0 * 1e6, sys.D_g[2] => 0.116e-6, sys.logσ[2] => 0.217,
        sys.N[3] => 0.300 * 1e6,  sys.D_g[3] => 1.800e-6, sys.logσ[3] => 0.380,
    )
    return System(equations(sys), t; name, defaults,
        systems = ModelingToolkit.get_systems(sys))
end

"""
    FreeTroposphereAerosol(; name=:FreeTroposphereAerosol)

Free troposphere aerosol size distribution from Table 8.3 of Seinfeld and Pandis (2006).

Three lognormal modes with parameters:
- Mode 1: (N=129 cm⁻³, Dg=0.007 μm, log σ=0.645)
- Mode 2: (N=59.7 cm⁻³, Dg=0.250 μm, log σ=0.253)
- Mode 3: (N=63.5 cm⁻³, Dg=0.520 μm, log σ=0.425)
"""
@component function FreeTroposphereAerosol(; name = :FreeTroposphereAerosol)
    sys = AerosolDistribution(3; name)
    defaults = Dict(
        sys.N[1] => 129.0 * 1e6, sys.D_g[1] => 0.007e-6, sys.logσ[1] => 0.645,
        sys.N[2] => 59.7 * 1e6,  sys.D_g[2] => 0.250e-6, sys.logσ[2] => 0.253,
        sys.N[3] => 63.5 * 1e6,  sys.D_g[3] => 0.520e-6, sys.logσ[3] => 0.425,
    )
    return System(equations(sys), t; name, defaults,
        systems = ModelingToolkit.get_systems(sys))
end

"""
    PolarAerosol(; name=:PolarAerosol)

Polar aerosol size distribution from Table 8.3 of Seinfeld and Pandis (2006).

Three lognormal modes with parameters:
- Mode 1: (N=21.7 cm⁻³, Dg=0.138 μm, log σ=0.164)
- Mode 2: (N=0.186 cm⁻³, Dg=0.750 μm, log σ=0.521)
- Mode 3: (N=3.04×10⁻⁴ cm⁻³, Dg=8.600 μm, log σ=0.420)
"""
@component function PolarAerosol(; name = :PolarAerosol)
    sys = AerosolDistribution(3; name)
    defaults = Dict(
        sys.N[1] => 21.7 * 1e6,    sys.D_g[1] => 0.138e-6, sys.logσ[1] => 0.164,
        sys.N[2] => 0.186 * 1e6,   sys.D_g[2] => 0.750e-6, sys.logσ[2] => 0.521,
        sys.N[3] => 3.04e-4 * 1e6, sys.D_g[3] => 8.600e-6, sys.logσ[3] => 0.420,
    )
    return System(equations(sys), t; name, defaults,
        systems = ModelingToolkit.get_systems(sys))
end

"""
    DesertAerosol(; name=:DesertAerosol)

Desert aerosol size distribution from Table 8.3 of Seinfeld and Pandis (2006).

Three lognormal modes with parameters:
- Mode 1: (N=726 cm⁻³, Dg=0.002 μm, log σ=0.247)
- Mode 2: (N=114 cm⁻³, Dg=0.038 μm, log σ=0.770)
- Mode 3: (N=0.178 cm⁻³, Dg=21.60 μm, log σ=0.438)
"""
@component function DesertAerosol(; name = :DesertAerosol)
    sys = AerosolDistribution(3; name)
    defaults = Dict(
        sys.N[1] => 726.0 * 1e6,  sys.D_g[1] => 0.002e-6, sys.logσ[1] => 0.247,
        sys.N[2] => 114.0 * 1e6,  sys.D_g[2] => 0.038e-6, sys.logσ[2] => 0.770,
        sys.N[3] => 0.178 * 1e6,  sys.D_g[3] => 21.60e-6, sys.logσ[3] => 0.438,
    )
    return System(equations(sys), t; name, defaults,
        systems = ModelingToolkit.get_systems(sys))
end
