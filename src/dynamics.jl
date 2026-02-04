export AerosolDynamics, DiameterGrowthRate, BrownianCoagulationCoefficient,
       MonodisperseCoagulation, DiscreteCoagulation

"""
    DiameterGrowthRate(; name=:DiameterGrowthRate)

Compute the diameter growth rate due to condensation in the continuum regime.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 13, Eq. 13.11-13.13.

The diameter growth rate is given by:

``I_D = \\frac{dD_p}{dt} = \\frac{4 D_i M_i}{R T D_p \\rho_p} (p_i - p_{eq,i})``

In the continuum regime with constant supersaturation (A = constant), this simplifies to:

``I_D = \\frac{A}{D_p}``

where ``A = \\frac{4 D_i M_i}{R T \\rho_p} (p_i - p_{eq,i})``.
"""
@component function DiameterGrowthRate(; name = :DiameterGrowthRate)
    @constants begin
        R_gas = 8.314, [description = "Universal gas constant", unit = u"J/mol/K"]
    end

    @parameters begin
        T = 298.0, [description = "Temperature", unit = u"K"]
        D_diff = 1e-5,
        [description = "Diffusion coefficient of condensing species", unit = u"m^2/s"]
        M_i = 0.1,
        [description = "Molecular weight of condensing species", unit = u"kg/mol"]
        ρ_p = 1000.0, [description = "Particle density", unit = u"kg/m^3"]
        Δp = 1e-4, [description = "Vapor pressure difference (p_i - p_eq)", unit = u"Pa"]
    end

    @variables begin
        D_p(t), [description = "Particle diameter", unit = u"m"]
        I_D(t), [description = "Diameter growth rate", unit = u"m/s"]
        A(t), [description = "Growth parameter", unit = u"m^2/s"]
    end

    eqs = [
        # Eq. 13.12 - Growth parameter (simplified form)
        A ~ 4 * D_diff * M_i * Δp / (R_gas * T * ρ_p),
        # Eq. 13.13 - Diameter growth rate
        I_D ~ A / D_p,
        # Eq. 13.20 - ODE for particle diameter evolution
        D(D_p) ~ I_D
    ]

    return System(eqs, t; name)
end

"""
    BrownianCoagulationCoefficient(; name=:BrownianCoagulationCoefficient)

Compute the Brownian coagulation coefficient using the Fuchs form, which is valid
across the continuum, transition, and free-molecular regimes.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 13, Table 13.1.

The Fuchs form of the coagulation coefficient is:

``K_{12} = 2\\pi(D_1 + D_2)(D_{p1} + D_{p2}) \\left( \\frac{D_{p1} + D_{p2}}{D_{p1} + D_{p2} + 2\\sqrt{g_1^2 + g_2^2}} + \\frac{8(D_1 + D_2)}{\\sqrt{\\bar{c}_1^2 + \\bar{c}_2^2}(D_{p1} + D_{p2})} \\right)^{-1}``

where:

  - ``D_i = \\frac{k T C_c}{3\\pi\\mu D_{pi}}`` is the Brownian diffusion coefficient
  - ``\\bar{c}_i = \\sqrt{\\frac{8 k T}{\\pi m_i}}`` is the mean thermal velocity
  - ``\\ell_i = \\frac{8 D_i}{\\pi \\bar{c}_i}`` is the particle mean free path
  - ``g_i = \\frac{1}{3 D_{pi} \\ell_i}\\left[(D_{pi} + \\ell_i)^3 - (D_{pi}^2 + \\ell_i^2)^{3/2}\\right] - D_{pi}``
"""
@component function BrownianCoagulationCoefficient(; name = :BrownianCoagulationCoefficient)
    @constants begin
        k_B = 1.380649e-23, [description = "Boltzmann constant", unit = u"J/K"]
        π_val = 3.14159265358979, [description = "Pi (dimensionless)", unit = u"1"]
    end

    @parameters begin
        T = 298.0, [description = "Temperature", unit = u"K"]
        μ = 1.83e-5, [description = "Dynamic viscosity of air", unit = u"Pa*s"]
        ρ_p = 1000.0, [description = "Particle density", unit = u"kg/m^3"]
        λ_air = 6.86e-8, [description = "Mean free path of air", unit = u"m"]
        D_p1 = 1e-7, [description = "Diameter of particle 1", unit = u"m"]
        D_p2 = 1e-7, [description = "Diameter of particle 2", unit = u"m"]
    end

    @variables begin
        # Particle masses (Eq. from volume and density)
        m_1(t), [description = "Mass of particle 1", unit = u"kg"]
        m_2(t), [description = "Mass of particle 2", unit = u"kg"]
        # Slip correction factors (Cunningham correction)
        C_c1(t),
        [description = "Slip correction factor for particle 1 (dimensionless)", unit = u"1"]
        C_c2(t),
        [description = "Slip correction factor for particle 2 (dimensionless)", unit = u"1"]
        Kn_1(t),
        [description = "Knudsen number for particle 1 (dimensionless)", unit = u"1"]
        Kn_2(t),
        [description = "Knudsen number for particle 2 (dimensionless)", unit = u"1"]
        # Brownian diffusion coefficients (Eq. 13.47)
        D_1(t), [description = "Diffusion coefficient of particle 1", unit = u"m^2/s"]
        D_2(t), [description = "Diffusion coefficient of particle 2", unit = u"m^2/s"]
        # Mean thermal velocities (Table 13.1)
        c_bar_1(t), [description = "Mean thermal velocity of particle 1", unit = u"m/s"]
        c_bar_2(t), [description = "Mean thermal velocity of particle 2", unit = u"m/s"]
        # Particle mean free paths (Table 13.1)
        ℓ_1(t), [description = "Mean free path of particle 1", unit = u"m"]
        ℓ_2(t), [description = "Mean free path of particle 2", unit = u"m"]
        # Fuchs interpolation parameters (Table 13.1)
        g_1(t), [description = "Fuchs interpolation parameter for particle 1", unit = u"m"]
        g_2(t), [description = "Fuchs interpolation parameter for particle 2", unit = u"m"]
        # Coagulation coefficient (Table 13.1)
        K_12(t), [description = "Brownian coagulation coefficient", unit = u"m^3/s"]
    end

    eqs = [
        # Particle masses from diameter and density
        m_1 ~ ρ_p * π_val / 6 * D_p1^3,
        m_2 ~ ρ_p * π_val / 6 * D_p2^3,

        # Knudsen numbers (Kn = 2λ/D_p)
        Kn_1 ~ 2 * λ_air / D_p1,
        Kn_2 ~ 2 * λ_air / D_p2,

        # Cunningham slip correction factor (standard form)
        C_c1 ~ 1 + Kn_1 * (1.257 + 0.4 * exp(-1.1 / Kn_1)),
        C_c2 ~ 1 + Kn_2 * (1.257 + 0.4 * exp(-1.1 / Kn_2)),

        # Eq. 13.47 - Brownian diffusion coefficient with slip correction (Table 13.1)
        D_1 ~ k_B * T * C_c1 / (3 * π_val * μ * D_p1),
        D_2 ~ k_B * T * C_c2 / (3 * π_val * μ * D_p2),

        # Table 13.1 - Mean thermal velocities
        c_bar_1 ~ sqrt(8 * k_B * T / (π_val * m_1)),
        c_bar_2 ~ sqrt(8 * k_B * T / (π_val * m_2)),

        # Table 13.1 - Particle mean free paths
        ℓ_1 ~ 8 * D_1 / (π_val * c_bar_1),
        ℓ_2 ~ 8 * D_2 / (π_val * c_bar_2),

        # Table 13.1 - Fuchs interpolation parameters
        g_1 ~ (1 / (3 * D_p1 * ℓ_1)) * ((D_p1 + ℓ_1)^3 - (D_p1^2 + ℓ_1^2)^1.5) - D_p1,
        g_2 ~ (1 / (3 * D_p2 * ℓ_2)) * ((D_p2 + ℓ_2)^3 - (D_p2^2 + ℓ_2^2)^1.5) - D_p2,

        # Table 13.1 - Fuchs form of coagulation coefficient
        K_12 ~
        2 * π_val * (D_1 + D_2) * (D_p1 + D_p2) / (
            (D_p1 + D_p2) / (D_p1 + D_p2 + 2 * sqrt(g_1^2 + g_2^2)) +
            8 * (D_1 + D_2) / (sqrt(c_bar_1^2 + c_bar_2^2) * (D_p1 + D_p2))
        )
    ]

    return System(eqs, t; name)
end

"""
    MonodisperseCoagulation(; name=:MonodisperseCoagulation)

Model the evolution of total aerosol number concentration for a monodisperse
aerosol undergoing Brownian coagulation with a constant coagulation coefficient.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 13, Eq. 13.65-13.67.

The governing equation is:

``\\frac{dN}{dt} = -\\frac{1}{2} K N^2(t)``

with analytical solution:

``N(t) = \\frac{N_0}{1 + t/\\tau_c}``

where ``\\tau_c = \\frac{2}{K N_0}`` is the characteristic coagulation time.
"""
@component function MonodisperseCoagulation(; name = :MonodisperseCoagulation)
    @parameters begin
        K = 1e-15, [description = "Coagulation coefficient", unit = u"m^3/s"]
        N_0 = 1e12, [description = "Initial number concentration", unit = u"m^-3"]
    end

    @variables begin
        N(t), [description = "Total number concentration", unit = u"m^-3"]
        τ_c(t), [description = "Characteristic coagulation time", unit = u"s"]
    end

    eqs = [
        # Eq. 13.67 - Characteristic coagulation time
        τ_c ~ 2 / (K * N_0),
        # Eq. 13.65 - Rate of change of total number concentration
        D(N) ~ -0.5 * K * N^2
    ]

    return System(eqs, t; name)
end

"""
    DiscreteCoagulation(n_bins::Int; name=:DiscreteCoagulation)

Model the discrete coagulation equation for a monodisperse initial distribution
with constant coagulation coefficient.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 13, Eq. 13.59, 13.71.

The discrete coagulation equation for k-mers is:

``\\frac{dN_k}{dt} = \\frac{1}{2} \\sum_{j=1}^{k-1} K_{j,k-j} N_j N_{k-j} - N_k \\sum_{j=1}^{\\infty} K_{kj} N_j``

For a constant coagulation coefficient K and monodisperse initial condition
(all particles initially monomers), the analytical solution is:

``N_k(t) = \\frac{N_0 (t/\\tau_c)^{k-1}}{[1 + t/\\tau_c]^{k+1}}``

Arguments:

  - `n_bins::Int`: Number of size bins (k-mer sizes to track)
"""
@component function DiscreteCoagulation(n_bins::Int; name = :DiscreteCoagulation)
    @parameters begin
        K = 1e-15, [description = "Coagulation coefficient", unit = u"m^3/s"]
        N_0 = 1e12, [description = "Initial total number concentration", unit = u"m^-3"]
    end

    @variables begin
        (N(t))[1:n_bins], [description = "Number concentration of k-mers", unit = u"m^-3"]
        τ_c(t), [description = "Characteristic coagulation time", unit = u"s"]
        N_total(t), [description = "Total number concentration", unit = u"m^-3"]
    end

    # Build the discrete coagulation equations (Eq. 13.59)
    # For constant K, this simplifies to:
    # dN_k/dt = (K/2) * sum_{j=1}^{k-1} N_j * N_{k-j} - K * N_k * sum_{j=1}^{n_bins} N_j

    coag_eqs = []
    for k in 1:n_bins
        # Formation term: particles of size j and k-j collide to form size k
        formation = k > 1 ? sum(N[j] * N[k - j] for j in 1:(k - 1)) : 0.0

        # Loss term: particle of size k collides with any other particle
        loss = N[k] * N_total

        push!(coag_eqs, D(N[k]) ~ 0.5 * K * formation - K * loss)
    end

    eqs = [
        # Eq. 13.67 - Characteristic coagulation time
        τ_c ~ 2 / (K * N_0),
        # Total number concentration
        N_total ~ sum(N[k] for k in 1:n_bins),
        # Coagulation equations for each bin
        coag_eqs...
    ]

    return System(eqs, t; name)
end

"""
    AerosolDynamics(; name=:AerosolDynamics)

A combined system for modeling aerosol population dynamics including condensation
growth and coagulation processes.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006) *Atmospheric Chemistry and Physics*,
2nd Edition, Chapter 13.

This system provides:

  - Condensation growth rate calculation
  - Brownian coagulation coefficient (Fuchs form)
  - Monodisperse coagulation dynamics

The General Dynamic Equation (Eq. 13.86) combines these processes:

``\\frac{\\partial n(v,t)}{\\partial t} = \\left(\\frac{\\partial n}{\\partial t}\\right)_{coag} + \\left(\\frac{\\partial n}{\\partial t}\\right)_{cond} + J_0\\delta(v-v_0) + S(v) - R(v)``
"""
@component function AerosolDynamics(; name = :AerosolDynamics)
    # Create subsystems
    growth = DiameterGrowthRate(; name = :growth)
    coag_coeff = BrownianCoagulationCoefficient(; name = :coag_coeff)
    coag = MonodisperseCoagulation(; name = :coag)

    # Connect subsystems through shared parameters
    @parameters begin
        T = 298.0, [description = "Temperature", unit = u"K"]
        ρ_p = 1000.0, [description = "Particle density", unit = u"kg/m^3"]
    end

    eqs = [
        # Share temperature between subsystems
        growth.T ~ T,
        coag_coeff.T ~ T,
        # Share particle density
        growth.ρ_p ~ ρ_p,
        coag_coeff.ρ_p ~ ρ_p,
        # Use computed coagulation coefficient for dynamics
        coag.K ~ coag_coeff.K_12
    ]

    return System(eqs, t; systems = [growth, coag_coeff, coag], name)
end
