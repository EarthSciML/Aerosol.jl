export SingleParticleDynamics, MeanFreePath, SlipCorrection, SettlingVelocity,
       BrownianDiffusion, ParticleMobility, ElectricalMobility, StokesNumber,
       AerodynamicDiameter

"""
    MeanFreePath(; name=:MeanFreePath)

Component for calculating the mean free path of air molecules.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 9, Eq. 9.6.

The mean free path is calculated from:
```math
\\lambda = \\frac{2\\mu}{p\\sqrt{8M/(\\pi R T)}}
```

where μ is dynamic viscosity, p is pressure, M is molecular weight, R is the gas constant,
and T is temperature.
"""
function MeanFreePath(; name=:MeanFreePath)
    @constants begin
        R = 8.314, [unit = u"J/(mol*K)", description = "Molar gas constant"]
        M_air = 0.02897, [unit = u"kg/mol", description = "Molecular weight of air"]
        π_val = π, [unit = u"1", description = "Pi (dimensionless)"]
    end

    @parameters begin
        T = 298.15, [unit = u"K", description = "Temperature"]
        P = 101325.0, [unit = u"Pa", description = "Pressure"]
        μ = 1.8e-5, [unit = u"kg/(m*s)", description = "Dynamic viscosity of air"]
    end

    @variables begin
        λ(t), [unit = u"m", description = "Mean free path of air"]
    end

    eqs = [
        # Eq. 9.6: λ = 2μ / [p * sqrt(8M/(πRT))]
        λ ~ 2 * μ / (P * sqrt(8 * M_air / (π_val * R * T)))
    ]

    System(eqs, t; name=name)
end

"""
    SlipCorrection(; name=:SlipCorrection)

Component for calculating the Cunningham slip correction factor.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 9, Eq. 9.34.

The slip correction factor accounts for noncontinuum effects when the particle size
approaches the mean free path of gas molecules:
```math
C_c = 1 + \\frac{2\\lambda}{D_p}\\left[A_1 + A_2 \\exp\\left(-\\frac{A_3 D_p}{2\\lambda}\\right)\\right]
```

where the constants A₁ = 1.257, A₂ = 0.4, A₃ = 1.1 are from Allen & Raabe (1982).
"""
function SlipCorrection(; name=:SlipCorrection)
    @constants begin
        A1 = 1.257, [unit = u"1", description = "Slip correction constant A1 (dimensionless)"]
        A2 = 0.4, [unit = u"1", description = "Slip correction constant A2 (dimensionless)"]
        A3 = 1.1, [unit = u"1", description = "Slip correction constant A3 (dimensionless)"]
    end

    @parameters begin
        D_p = 1e-6, [unit = u"m", description = "Particle diameter"]
        λ = 6.51e-8, [unit = u"m", description = "Mean free path of air"]
    end

    @variables begin
        C_c(t), [unit = u"1", description = "Cunningham slip correction factor (dimensionless)"]
        Kn(t), [unit = u"1", description = "Knudsen number (dimensionless)"]
    end

    eqs = [
        # Eq. 9.1: Kn = 2λ/D_p
        Kn ~ 2 * λ / D_p,
        # Eq. 9.34: C_c = 1 + Kn * [A1 + A2 * exp(-A3/Kn)]
        C_c ~ 1 + Kn * (A1 + A2 * exp(-A3 / Kn))
    ]

    System(eqs, t; name=name)
end

"""
    SettlingVelocity(; name=:SettlingVelocity)

Component for calculating particle terminal settling velocity and relaxation time.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 9, Eqs. 9.38, 9.42.

The terminal settling velocity is given by:
```math
v_t = \\frac{D_p^2 \\rho_p g C_c}{18\\mu}
```

The relaxation time (time constant for velocity adjustment) is:
```math
\\tau = \\frac{\\rho_p D_p^2 C_c}{18\\mu}
```
"""
function SettlingVelocity(; name=:SettlingVelocity)
    @constants begin
        g = 9.807, [unit = u"m/s^2", description = "Gravitational acceleration"]
    end

    @parameters begin
        D_p = 1e-6, [unit = u"m", description = "Particle diameter"]
        ρ_p = 1000.0, [unit = u"kg/m^3", description = "Particle density"]
        μ = 1.8e-5, [unit = u"kg/(m*s)", description = "Dynamic viscosity of air"]
        C_c = 1.0, [unit = u"1", description = "Slip correction factor (dimensionless)"]
    end

    @variables begin
        v_t(t), [unit = u"m/s", description = "Terminal settling velocity"]
        τ(t), [unit = u"s", description = "Relaxation time"]
    end

    eqs = [
        # Eq. 9.42: v_t = D_p² ρ_p g C_c / (18 μ)
        v_t ~ (D_p^2 * ρ_p * g * C_c) / (18 * μ),
        # Eq. 9.38: τ = ρ_p D_p² C_c / (18 μ)
        τ ~ (ρ_p * D_p^2 * C_c) / (18 * μ)
    ]

    System(eqs, t; name=name)
end

"""
    BrownianDiffusion(; name=:BrownianDiffusion)

Component for calculating Brownian diffusion properties of aerosol particles.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 9, Eqs. 9.73, 9.78, 9.87, 9.89.

The Stokes-Einstein-Sutherland relation gives the diffusion coefficient:
```math
D = \\frac{k T C_c}{3 \\pi \\mu D_p}
```

The particle mobility is:
```math
B = \\frac{C_c}{3 \\pi \\mu D_p}
```

The mean thermal speed of the particle is:
```math
\\bar{c}_p = \\sqrt{\\frac{8 k T}{\\pi m_p}}
```

The particle mean free path is derived from D = ½ c̄_p λ_p (Eq. 9.88):
```math
\\lambda_p = \\frac{2D}{\\bar{c}_p}
```
"""
function BrownianDiffusion(; name=:BrownianDiffusion)
    @constants begin
        k_B = 1.381e-23, [unit = u"J/K", description = "Boltzmann constant"]
        π_val = π, [unit = u"1", description = "Pi (dimensionless)"]
    end

    @parameters begin
        D_p = 1e-6, [unit = u"m", description = "Particle diameter"]
        ρ_p = 1000.0, [unit = u"kg/m^3", description = "Particle density"]
        T = 298.15, [unit = u"K", description = "Temperature"]
        μ = 1.8e-5, [unit = u"kg/(m*s)", description = "Dynamic viscosity of air"]
        C_c = 1.0, [unit = u"1", description = "Slip correction factor (dimensionless)"]
    end

    @variables begin
        D_B(t), [unit = u"m^2/s", description = "Brownian diffusion coefficient"]
        B(t), [unit = u"s/kg", description = "Particle mobility"]
        c_p(t), [unit = u"m/s", description = "Mean thermal speed of particle"]
        λ_p(t), [unit = u"m", description = "Particle mean free path"]
        m_p(t), [unit = u"kg", description = "Particle mass"]
    end

    eqs = [
        # Particle mass from density and volume
        m_p ~ ρ_p * (π_val / 6) * D_p^3,
        # Eq. 9.73: D = k T C_c / (3 π μ D_p)
        D_B ~ (k_B * T * C_c) / (3 * π_val * μ * D_p),
        # Eq. 9.78: B = C_c / (3 π μ D_p)
        B ~ C_c / (3 * π_val * μ * D_p),
        # Eq. 9.87: c̄_p = sqrt(8 k T / (π m_p))
        c_p ~ sqrt(8 * k_B * T / (π_val * m_p)),
        # Eq. 9.88: λ_p = 2D / c̄_p (from D = ½ c̄_p λ_p)
        λ_p ~ 2 * D_B / c_p
    ]

    System(eqs, t; name=name)
end

"""
    ParticleMobility(; name=:ParticleMobility)

Component for calculating particle mechanical mobility.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 9, Eq. 9.78.

The particle mobility relates external force to drift velocity:
```math
B = \\frac{C_c}{3 \\pi \\mu D_p}
```

With the Einstein relation connecting mobility to diffusion:
```math
D = B k T
```
"""
function ParticleMobility(; name=:ParticleMobility)
    @constants begin
        π_val = π, [unit = u"1", description = "Pi (dimensionless)"]
    end

    @parameters begin
        D_p = 1e-6, [unit = u"m", description = "Particle diameter"]
        μ = 1.8e-5, [unit = u"kg/(m*s)", description = "Dynamic viscosity of air"]
        C_c = 1.0, [unit = u"1", description = "Slip correction factor (dimensionless)"]
    end

    @variables begin
        B(t), [unit = u"s/kg", description = "Particle mobility"]
    end

    eqs = [
        # Eq. 9.78: B = C_c / (3 π μ D_p)
        B ~ C_c / (3 * π_val * μ * D_p)
    ]

    System(eqs, t; name=name)
end

"""
    ElectricalMobility(; name=:ElectricalMobility)

Component for calculating electrical mobility and migration velocity.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 9, Eqs. 9.49, 9.50, 9.51.

The electrical mobility relates particle velocity to electric field strength:
```math
B_e = \\frac{q C_c}{3 \\pi \\mu D_p}
```

The migration velocity in an electric field E is:
```math
v_e = B_e E
```
"""
function ElectricalMobility(; name=:ElectricalMobility)
    @constants begin
        π_val = π, [unit = u"1", description = "Pi (dimensionless)"]
        e = 1.602e-19, [unit = u"C", description = "Elementary charge"]
    end

    @parameters begin
        D_p = 1e-6, [unit = u"m", description = "Particle diameter"]
        μ = 1.8e-5, [unit = u"kg/(m*s)", description = "Dynamic viscosity of air"]
        C_c = 1.0, [unit = u"1", description = "Slip correction factor (dimensionless)"]
        n_charges = 1.0, [unit = u"1", description = "Number of elementary charges (dimensionless)"]
        E = 1000.0, [unit = u"V/m", description = "Electric field strength"]
    end

    @variables begin
        q(t), [unit = u"C", description = "Particle charge"]
        B_e(t), [unit = u"m^2/(V*s)", description = "Electrical mobility"]
        v_e(t), [unit = u"m/s", description = "Electrical migration velocity"]
    end

    eqs = [
        # Particle charge
        q ~ n_charges * e,
        # Eq. 9.50: B_e = q C_c / (3 π μ D_p)
        B_e ~ (q * C_c) / (3 * π_val * μ * D_p),
        # Eq. 9.51: v_e = B_e E
        v_e ~ B_e * E
    ]

    System(eqs, t; name=name)
end

"""
    StokesNumber(; name=:StokesNumber)

Component for calculating Stokes number and stop distance.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 9, Eqs. 9.98, 9.101.

The stop distance is the distance a particle travels before coming to rest:
```math
s_p = \\tau U
```

The Stokes number characterizes particle inertia relative to flow:
```math
St = \\frac{\\tau u_0}{L} = \\frac{D_p^2 \\rho_p C_c u_0}{18 \\mu L}
```

When St ≫ 1, particles cannot follow fluid streamlines (inertial impaction).
When St ≪ 1, particles follow fluid motion closely.
"""
function StokesNumber(; name=:StokesNumber)
    @parameters begin
        D_p = 1e-6, [unit = u"m", description = "Particle diameter"]
        ρ_p = 1000.0, [unit = u"kg/m^3", description = "Particle density"]
        μ = 1.8e-5, [unit = u"kg/(m*s)", description = "Dynamic viscosity of air"]
        C_c = 1.0, [unit = u"1", description = "Slip correction factor (dimensionless)"]
        U = 1.0, [unit = u"m/s", description = "Initial particle velocity"]
        L = 0.01, [unit = u"m", description = "Characteristic length scale"]
    end

    @variables begin
        τ(t), [unit = u"s", description = "Relaxation time"]
        s_p(t), [unit = u"m", description = "Stop distance"]
        St(t), [unit = u"1", description = "Stokes number (dimensionless)"]
    end

    eqs = [
        # Eq. 9.38: τ = ρ_p D_p² C_c / (18 μ)
        τ ~ (ρ_p * D_p^2 * C_c) / (18 * μ),
        # Eq. 9.98: s_p = τ U
        s_p ~ τ * U,
        # Eq. 9.101: St = τ u_0 / L
        St ~ τ * U / L
    ]

    System(eqs, t; name=name)
end

"""
    AerodynamicDiameter(; name=:AerodynamicDiameter)

Component for calculating equivalent particle diameters.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 9, Eqs. 9.102, 9.106, 9.108, 9.110, 9.114.

The volume equivalent diameter is the diameter of a sphere with the same volume:
```math
D_{ve} = \\left(\\frac{6 V_p}{\\pi}\\right)^{1/3}
```

The aerodynamic diameter is defined for a sphere of unit density (1000 kg/m³)
with the same settling velocity:
```math
D_{ca} = D_p \\sqrt{\\frac{\\rho_p}{\\rho_p^\\circ}} \\sqrt{\\frac{C_c(D_p)}{C_c(D_{ca})}}
```

In the continuum regime (large particles), this simplifies to:
```math
D_{ca} \\approx D_p \\sqrt{\\frac{\\rho_p}{\\rho_p^\\circ}}
```

The electrical mobility equivalent diameter relates to mobility measurement:
```math
D_{em} = D_{ve} \\chi \\frac{C_c(D_{em})}{C_c(D_{ve})}
```
"""
function AerodynamicDiameter(; name=:AerodynamicDiameter)
    @constants begin
        ρ_p0 = 1000.0, [unit = u"kg/m^3", description = "Reference density (1 g/cm³)"]
    end

    @parameters begin
        D_p = 1e-6, [unit = u"m", description = "Particle diameter"]
        ρ_p = 1500.0, [unit = u"kg/m^3", description = "Particle density"]
        C_c_Dp = 1.0, [unit = u"1", description = "Slip correction at D_p (dimensionless)"]
        C_c_Da = 1.0, [unit = u"1", description = "Slip correction at D_ca (dimensionless)"]
        χ = 1.0, [unit = u"1", description = "Dynamic shape factor (dimensionless)"]
    end

    @variables begin
        D_ca(t), [unit = u"m", description = "Classical aerodynamic diameter"]
        D_ca_continuum(t), [unit = u"m", description = "Aerodynamic diameter (continuum limit)"]
    end

    eqs = [
        # Eq. 9.110: D_ca = D_p * sqrt(ρ_p/ρ_p0) * sqrt(C_c(D_p)/C_c(D_ca))
        D_ca ~ D_p * sqrt(ρ_p / ρ_p0) * sqrt(C_c_Dp / C_c_Da),
        # Eq. 9.112: Continuum limit (C_c ≈ 1)
        D_ca_continuum ~ D_p * sqrt(ρ_p / ρ_p0)
    ]

    System(eqs, t; name=name)
end

"""
    SingleParticleDynamics(; name=:SingleParticleDynamics)

Comprehensive component combining all single particle dynamics equations.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 9: Dynamics of Single Aerosol Particles.

This component calculates key aerosol particle properties:
- Mean free path of air (Eq. 9.6)
- Knudsen number (Eq. 9.1)
- Cunningham slip correction factor (Eq. 9.34)
- Terminal settling velocity (Eq. 9.42)
- Relaxation time (Eq. 9.38)
- Brownian diffusion coefficient (Eq. 9.73)
- Particle mobility (Eq. 9.78)
- Mean thermal speed (Eq. 9.87)
- Particle mean free path (Eq. 9.89)
- Aerodynamic diameter (Eq. 9.110)
"""
function SingleParticleDynamics(; name=:SingleParticleDynamics)
    @constants begin
        k_B = 1.381e-23, [unit = u"J/K", description = "Boltzmann constant"]
        R = 8.314, [unit = u"J/(mol*K)", description = "Molar gas constant"]
        g = 9.807, [unit = u"m/s^2", description = "Gravitational acceleration"]
        M_air = 0.02897, [unit = u"kg/mol", description = "Molecular weight of air"]
        π_val = π, [unit = u"1", description = "Pi (dimensionless)"]
        ρ_p0 = 1000.0, [unit = u"kg/m^3", description = "Reference density (1 g/cm³)"]
        # Slip correction constants (Allen & Raabe 1982)
        A1 = 1.257, [unit = u"1", description = "Slip correction constant A1 (dimensionless)"]
        A2 = 0.4, [unit = u"1", description = "Slip correction constant A2 (dimensionless)"]
        A3 = 1.1, [unit = u"1", description = "Slip correction constant A3 (dimensionless)"]
    end

    @parameters begin
        D_p = 1e-6, [unit = u"m", description = "Particle diameter"]
        ρ_p = 1000.0, [unit = u"kg/m^3", description = "Particle density"]
        T = 298.15, [unit = u"K", description = "Temperature"]
        P = 101325.0, [unit = u"Pa", description = "Pressure"]
        μ = 1.8e-5, [unit = u"kg/(m*s)", description = "Dynamic viscosity of air"]
    end

    @variables begin
        # Mean free path and Knudsen number
        λ(t), [unit = u"m", description = "Mean free path of air"]
        Kn(t), [unit = u"1", description = "Knudsen number (dimensionless)"]
        # Slip correction
        C_c(t), [unit = u"1", description = "Cunningham slip correction factor (dimensionless)"]
        # Settling
        v_t(t), [unit = u"m/s", description = "Terminal settling velocity"]
        τ(t), [unit = u"s", description = "Relaxation time"]
        # Brownian motion
        D_B(t), [unit = u"m^2/s", description = "Brownian diffusion coefficient"]
        B(t), [unit = u"s/kg", description = "Particle mobility"]
        m_p(t), [unit = u"kg", description = "Particle mass"]
        c_p(t), [unit = u"m/s", description = "Mean thermal speed of particle"]
        λ_p(t), [unit = u"m", description = "Particle mean free path"]
        # Aerodynamic diameter (continuum approximation)
        D_ca(t), [unit = u"m", description = "Aerodynamic diameter (continuum limit)"]
    end

    eqs = [
        # Eq. 9.6: Mean free path of air
        λ ~ 2 * μ / (P * sqrt(8 * M_air / (π_val * R * T))),
        # Eq. 9.1: Knudsen number
        Kn ~ 2 * λ / D_p,
        # Eq. 9.34: Cunningham slip correction factor
        C_c ~ 1 + Kn * (A1 + A2 * exp(-A3 / Kn)),
        # Particle mass
        m_p ~ ρ_p * (π_val / 6) * D_p^3,
        # Eq. 9.38: Relaxation time
        τ ~ (ρ_p * D_p^2 * C_c) / (18 * μ),
        # Eq. 9.42: Terminal settling velocity
        v_t ~ τ * g,
        # Eq. 9.73: Brownian diffusion coefficient
        D_B ~ (k_B * T * C_c) / (3 * π_val * μ * D_p),
        # Eq. 9.78: Particle mobility
        B ~ C_c / (3 * π_val * μ * D_p),
        # Eq. 9.87: Mean thermal speed
        c_p ~ sqrt(8 * k_B * T / (π_val * m_p)),
        # Eq. 9.88: Particle mean free path (from D = ½ c̄_p λ_p)
        λ_p ~ 2 * D_B / c_p,
        # Eq. 9.112: Aerodynamic diameter (continuum limit)
        D_ca ~ D_p * sqrt(ρ_p / ρ_p0)
    ]

    System(eqs, t; name=name)
end
