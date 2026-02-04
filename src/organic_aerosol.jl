export ECTracerMethod, NoninteractingSOA, AbsorptivePartitioning, TwoProductSOA,
       LangmuirAdsorption, BETAdsorption, FHHAdsorption

"""
    ECTracerMethod(; name=:ECTracerMethod)

The elemental carbon (EC) tracer method for estimating secondary organic carbon (OC)
concentrations from measurements of total OC and EC.

From Seinfeld and Pandis (2006), Chapter 14, Section 14.3.2, Equations 14.1-14.4.

Primary OC is estimated as the sum of combustion-related OC (proportional to EC via
the `OC_EC_ratio`) and noncombustion OC (`OC_NC`). Secondary OC is the remainder.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 14. John Wiley & Sons, Inc.
"""
@component function ECTracerMethod(; name = :ECTracerMethod)
    @constants begin
        zero_conc = 0.0,
        [description = "Zero concentration for unit consistency", unit = u"kg/m^3"]
    end

    @parameters begin
        OC_EC_ratio = 1.7,
        [description = "OC-to-EC ratio for combustion sources (dimensionless)", unit = u"1"]
        OC_NC = 0.9e-9,
        [description = "Noncombustion primary OC concentration", unit = u"kg/m^3"]
    end

    @variables begin
        OC(t), [description = "Total organic carbon concentration", unit = u"kg/m^3"]
        EC(t), [description = "Elemental carbon concentration", unit = u"kg/m^3"]
        OC_primary(t),
        [description = "Primary organic carbon concentration", unit = u"kg/m^3"]
        OC_secondary(t),
        [description = "Secondary organic carbon concentration", unit = u"kg/m^3"]
    end

    eqs = [
        OC_primary ~ OC_EC_ratio * EC + OC_NC,               # Eq. 14.4 - Primary OC from EC tracer
        OC_secondary ~ max(zero_conc, OC - OC_primary)       # Eq. 14.3 - Secondary OC (non-negative)
    ]

    return System(eqs, t; name)
end

"""
    NoninteractingSOA(; name=:NoninteractingSOA)

Gas-particle partitioning model for a secondary organic aerosol (SOA) compound
that does not interact with preexisting aerosol components. The compound exists
in the aerosol phase as pure compound, so gas-aerosol equilibrium is reached
when the gas-phase concentration equals the saturation concentration.

From Seinfeld and Pandis (2006), Chapter 14, Section 14.5.2, Equations 14.6-14.18.

The aerosol concentration is zero until the reacted ROG exceeds the threshold
`ΔROG*_i = p°_i M_ROG / (a_i R T)`, after which the excess condenses.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 14. John Wiley & Sons, Inc.
"""
@component function NoninteractingSOA(; name = :NoninteractingSOA)
    @constants begin
        R_gas = 8.314, [description = "Universal gas constant", unit = u"Pa*m^3/(mol*K)"]
        zero_conc = 0.0,
        [description = "Zero concentration for unit consistency", unit = u"kg/m^3"]
        zero_dimless = 0.0,
        [description = "Zero dimensionless for unit consistency", unit = u"1"]
    end

    @parameters begin
        a_i = 0.05, [description = "Molar yield of product i (dimensionless)", unit = u"1"]
        p_i = 1.01325e-5, [description = "Vapor pressure of pure product i", unit = u"Pa"]
        M_i = 0.180, [description = "Molecular weight of product i", unit = u"kg/mol"]
        M_ROG = 0.150, [description = "Molecular weight of parent ROG", unit = u"kg/mol"]
        T = 298.0, [description = "Temperature", unit = u"K"]
    end

    @variables begin
        ΔROG(t), [description = "Reacted ROG concentration", unit = u"kg/m^3"]
        c_eq(t),
        [description = "Equilibrium gas-phase saturation concentration", unit = u"kg/m^3"]
        c_total(t), [description = "Total product concentration", unit = u"kg/m^3"]
        c_aer(t), [description = "Aerosol-phase concentration of product", unit = u"kg/m^3"]
        c_gas(t), [description = "Gas-phase concentration of product", unit = u"kg/m^3"]
        ΔROG_threshold(t),
        [description = "Threshold reacted ROG for SOA formation", unit = u"kg/m^3"]
        X_p(t),
        [
            description = "Mass fraction of product in aerosol phase (dimensionless)", unit = u"1"]
        Y(t), [description = "Aerosol mass yield (dimensionless)", unit = u"1"]
    end

    eqs = [
        c_eq ~ p_i * M_i / (R_gas * T),                                          # Eq. 14.7
        c_total ~ a_i * (M_i / M_ROG) * ΔROG,                                    # Eq. 14.11
        ΔROG_threshold ~ p_i * M_ROG / (a_i * R_gas * T),                         # Eq. 14.13
        c_aer ~ max(zero_conc, c_total - c_eq),                                   # Eq. 14.9/14.10
        c_gas ~ c_total - c_aer,                                                  # Eq. 14.6 rearranged
        X_p ~ ifelse(c_total > zero_conc, c_aer / c_total, zero_dimless),          # Eq. 14.14
        Y ~ ifelse(ΔROG > zero_conc, c_aer / ΔROG, zero_dimless)                 # Eq. 14.16
    ]

    return System(eqs, t; name)
end

"""
    AbsorptivePartitioning(; name=:AbsorptivePartitioning)

Gas-particle partitioning model for SOA formation by absorption into a preexisting
organic aerosol phase, assuming an ideal solution. There is no threshold for SOA
formation; secondary aerosol forms as soon as the ROG starts reacting.

From Seinfeld and Pandis (2006), Chapter 14, Section 14.5.2, Equations 14.19-14.26.

Uses the dilute approximation `c_{aer,i} << m_0` (Eq. 14.23), valid when no
single compound dominates the organic aerosol composition.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 14. John Wiley & Sons, Inc.
"""
@component function AbsorptivePartitioning(; name = :AbsorptivePartitioning)
    @constants begin
        R_gas = 8.314, [description = "Universal gas constant", unit = u"Pa*m^3/(mol*K)"]
        zero_conc = 0.0,
        [description = "Zero concentration for unit consistency", unit = u"kg/m^3"]
        zero_dimless = 0.0,
        [description = "Zero dimensionless for unit consistency", unit = u"1"]
    end

    @parameters begin
        a_i = 0.05, [description = "Molar yield of product i (dimensionless)", unit = u"1"]
        p_i = 1.01325e-5, [description = "Vapor pressure of pure product i", unit = u"Pa"]
        M_i = 0.180, [description = "Molecular weight of product i", unit = u"kg/mol"]
        M_ROG = 0.150, [description = "Molecular weight of parent ROG", unit = u"kg/mol"]
        M_0 = 0.200,
        [description = "Average molecular weight of preexisting organic aerosol",
            unit = u"kg/mol"]
        m_0 = 10.0e-9,
        [description = "Preexisting organic aerosol mass concentration", unit = u"kg/m^3"]
        T = 298.0, [description = "Temperature", unit = u"K"]
    end

    @variables begin
        ΔROG(t), [description = "Reacted ROG concentration", unit = u"kg/m^3"]
        c_aer(t),
        [description = "Aerosol-phase concentration of product i", unit = u"kg/m^3"]
        X_p(t),
        [
            description = "Mass fraction of product in aerosol phase (dimensionless)", unit = u"1"]
        Y(t), [description = "Aerosol mass yield (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 14.24 - Aerosol concentration with preexisting OA (dilute approximation)
        c_aer ~
        (a_i * R_gas * T / M_ROG) * (M_i * m_0 / (m_0 * R_gas * T + M_0 * p_i)) * ΔROG,
        # Eq. 14.25 - Mass fraction in aerosol phase
        X_p ~ m_0 * R_gas * T / (m_0 * R_gas * T + p_i * M_0),
        # Eq. 14.26 - Aerosol mass yield
        Y ~ ifelse(ΔROG > zero_conc, c_aer / ΔROG, zero_dimless)
    ]

    return System(eqs, t; name)
end

"""
    TwoProductSOA(; name=:TwoProductSOA)

The two-product model for secondary organic aerosol (SOA) formation, developed by
Odum et al. (1996). This model represents the complex mixture of oxidation products
from a reactive organic gas (ROG) as two lumped surrogate products that form an
ideal solution.

From Seinfeld and Pandis (2006), Chapter 14, Section 14.5.2, Equations 14.27-14.43.

The key equation (14.43) gives the overall aerosol mass yield:

    Y = c_aer * (a₁/(c_aer + c₁°) + a₂/(c_aer + c₂°))

where `a_i` are stoichiometric yields and `c_i°` are saturation concentrations.

Parameters from Table 14.12 are provided for multiple VOC precursors.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 14. John Wiley & Sons, Inc.

See also: Odum, J. R. et al. (1996). Gas/Particle Partitioning and Secondary Organic
Aerosol Yields. *Environ. Sci. Technol.*, 30, 2580-2585.
"""
@component function TwoProductSOA(; name = :TwoProductSOA)
    @parameters begin
        a_1 = 0.038,
        [
            description = "Stoichiometric mass yield of product 1 (dimensionless)", unit = u"1"]
        a_2 = 0.326,
        [
            description = "Stoichiometric mass yield of product 2 (dimensionless)", unit = u"1"]
        c_sat_1 = 5.8e-9,
        [description = "Saturation concentration of product 1", unit = u"kg/m^3"]
        c_sat_2 = 250.0e-9,
        [description = "Saturation concentration of product 2", unit = u"kg/m^3"]
    end

    @variables begin
        ΔROG(t), [description = "Reacted ROG concentration", unit = u"kg/m^3"]
        c_aer(t),
        [description = "Total organic aerosol concentration produced", unit = u"kg/m^3"]
        Y(t), [description = "Aerosol mass yield (dimensionless)", unit = u"1"]
        ΔROG_threshold(t),
        [description = "Threshold reacted ROG for binary solution SOA formation",
            unit = u"kg/m^3"]
    end

    eqs = [
        # Eq. 14.43 - Two-product Odum model (key equation)
        Y ~ c_aer * (a_1 / (c_aer + c_sat_1) + a_2 / (c_aer + c_sat_2)),
        # Eq. 14.39 - Threshold for binary solution formation
        ΔROG_threshold ~ 1 / (a_1 / c_sat_1 + a_2 / c_sat_2),
        # Mass balance: c_aer = Y * ΔROG
        c_aer ~ Y * ΔROG
    ]

    return System(eqs, t; name)
end

"""
    LangmuirAdsorption(; name=:LangmuirAdsorption)

Langmuir adsorption isotherm for the formation of a monolayer of adsorbed
molecules on particle surfaces. Assumes all adsorption sites are equivalent,
no horizontal interactions, and uniform heat of adsorption.

From Seinfeld and Pandis (2006), Chapter 14, Section 14.5.3, Equation 14.44.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 14. John Wiley & Sons, Inc.
"""
@component function LangmuirAdsorption(; name = :LangmuirAdsorption)
    @parameters begin
        V_m = 1.0, [description = "Gas volume for monolayer formation", unit = u"m^3"]
        b = 1.0, [description = "Adsorption constant", unit = u"Pa^-1"]
    end

    @variables begin
        p(t), [description = "Gas partial pressure", unit = u"Pa"]
        V(t), [description = "Volume of gas adsorbed", unit = u"m^3"]
    end

    eqs = [
        V ~ V_m * b * p / (1 + b * p),  # Eq. 14.44
    ]

    return System(eqs, t; name)
end

"""
    BETAdsorption(; name=:BETAdsorption)

BET (Brunauer-Emmett-Teller) adsorption isotherm for multilayer adsorption.
Extends the Langmuir model to include formation of additional molecular layers,
assuming each adsorbed molecule can serve as an adsorption site for the next layer.

From Seinfeld and Pandis (2006), Chapter 14, Section 14.5.3, Equation 14.45.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 14. John Wiley & Sons, Inc.
"""
@component function BETAdsorption(; name = :BETAdsorption)
    @parameters begin
        V_m = 1.0, [description = "Gas volume for monolayer formation", unit = u"m^3"]
        c_BET = 10.0,
        [
            description = "BET constant for the adsorbing surface (dimensionless)", unit = u"1"]
    end

    @variables begin
        S(t),
        [description = "Gas-phase saturation ratio p/p_sat (dimensionless)", unit = u"1"]
        V(t), [description = "Volume of gas adsorbed", unit = u"m^3"]
    end

    eqs = [
    # Eq. 14.45 - BET isotherm
        V ~ V_m * c_BET * S / ((1 - S) * (1 + (c_BET - 1) * S)),
    ]

    return System(eqs, t; name)
end

"""
    FHHAdsorption(; name=:FHHAdsorption)

FHH (Frenkel-Halsey-Hill) adsorption isotherm, appropriate for conditions
near saturation where the BET theory breaks down. Relates the logarithm of
the supersaturation ratio to the adsorbed volume.

From Seinfeld and Pandis (2006), Chapter 14, Section 14.5.3, Equation 14.46.

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 14. John Wiley & Sons, Inc.
"""
@component function FHHAdsorption(; name = :FHHAdsorption)
    @parameters begin
        V_m = 1.0, [description = "Gas volume for monolayer formation", unit = u"m^3"]
        A_FHH = 1.0, [description = "FHH constant A (dimensionless)", unit = u"1"]
        B_FHH = 1.0, [description = "FHH constant B (dimensionless)", unit = u"1"]
    end

    @variables begin
        S(t),
        [description = "Gas-phase saturation ratio p/p_sat (dimensionless)", unit = u"1"]
        V(t), [description = "Volume of gas adsorbed", unit = u"m^3"]
    end

    eqs = [
    # Eq. 14.46 - FHH isotherm: ln(p°/p) = A / (V/V_m)^B
    # Rearranged: V/V_m = (A / ln(1/S))^(1/B) since S = p/p°
        V ~ V_m * (A_FHH / log(1 / S))^(1 / B_FHH),
    ]

    return System(eqs, t; name)
end
