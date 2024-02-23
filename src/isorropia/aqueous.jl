"""
An ion with the given concentration `m` [mol/m³ air] and valence (charge) `z`.
"""
struct Ion <: Species
    m::Num
    z::Int
end

function Base.nameof(i::Ion)
    string(Symbolics.tosymbol(i.m, escape=false))
end

vars(i::Ion) = [i.m]
terms(i::Ion) = [i.m], [1]

#==
The activity coefficient of an ion is assumed to be one (Fountoukis and Nenes (2007), Section 3.3).
==#
activity(i::Ion, W) = i.m / W

"""
Generate the aqueous ions.
Each ion has an associated MTK variable named 
<name>_aq, where <name> is the name of the compound.
"""
function generate_ions(t)
    ion_names = [:NH4, :Na, :H, :Cl, :NO3, :SO4, :HNO3, :NH3, :HCl, :HSO4, :Ca, :K, :Mg, :OH]
    ion_valence = [1, 1, 1, -1, -1, -2, 0, 0, 0, -1, 2, 1, 2, -1]
    ions = Dict()
    for i ∈ eachindex(ion_names)
        varname = Symbol(ion_names[i], "_aq")
        s = "Aqueous $(ion_names[i])"
        var = only(@variables $varname(t) = 1e-8 [unit = u"mol/m_air^3", description = s])
        var = ParentScope(var)
        ions[ion_names[i]] = Ion(var, ion_valence[i])
    end
    ions
end

abstract type SaltLike <: Species end

"""
An aqueous salt comprised of a cation, an anion, and an activity parameter (q).
q values are given in Table 4 of Fountoukis and Nenes (2007).
"""
struct Salt <: SaltLike
    cation::Ion
    "Number of cations per molecule"
    ν_cation::Number
    anion::Ion
    "Number of anions per molecule"
    ν_anion::Number
    "Deliquescence relative humidity at 298.15K"
    drh
    "Enthalpy term (-18/1000R L_s m_s)"
    l_term
    "Kusik-Meissner Binary activity parameter"
    q::Number

    function Salt(cation::Ion, ν_cation::Number, anion::Ion, ν_anion::Number, drh, l_term, q::Number)
        if cation.z * ν_cation + anion.z * ν_anion ≠ 0
            if q != Inf # Special case for NH43HSO42, which doesn't balance.
                throw(ArgumentError("The charge of the cation and anion must sum to zero."))
            end
        end
        new(cation, ν_cation, anion, ν_anion, drh, l_term, q)
    end
end

function Base.nameof(s::SaltLike)
    c = replace(string(Symbolics.tosymbol(s.cation.m, escape=false)), "_aq" => "")
    a = replace(string(Symbolics.tosymbol(s.anion.m, escape=false)), "_aq" => "")
    "$(c)$(s.ν_cation > 1 ? s.ν_cation : "")$(a)$(s.ν_anion > 1 ? s.ν_anion : "")_aqs"
end

vars(s::SaltLike) = [s.cation.m, s.anion.m]
terms(s::SaltLike) = [s.cation.m, s.anion.m], [s.ν_cation, s.ν_anion]

"Generate Salts from Table 4, where `i` is a dictionary of ions."
function generate_salts(i) 
    Dict(
        :CaNO32 => Salt(i[:Ca], 1, i[:NO3], 2, 0.4906, 509.4, 0.93),
        :CaCl2 => Salt(i[:Ca], 1, i[:Cl], 2, 0.2830, 551.1, 2.4),
        :CaSO4 => CaSO4aqs(i[:Ca], 1, i[:SO4], 1, 0.9700, missing, NaN),
        :KHSO4 => KHSO4aqs(i[:K], 1, i[:HSO4], 1, 0.8600, missing, NaN),
        :K2SO4 => Salt(i[:K], 2, i[:SO4], 1, 0.9751, 35.6, -0.25),
        :KNO3 => Salt(i[:K], 1, i[:NO3], 1, 0.9248, missing, -2.33),
        :KCl => Salt(i[:K], 1, i[:Cl], 1, 0.8426, 158.9, 0.92),
        :MgSO4 => Salt(i[:Mg], 1, i[:SO4], 1, 0.8613, -714.5, 0.15),
        :MgNO32 => Salt(i[:Mg], 1, i[:NO3], 2, 0.5400, 230.2, 2.32),
        :MgCl2 => Salt(i[:Mg], 1, i[:Cl], 2, 0.3284, 42.23, 2.90),
        :NaCl => Salt(i[:Na], 1, i[:Cl], 1, 0.7528, 25.0, 2.23),
        :Na2SO4 => Salt(i[:Na], 2, i[:SO4], 1, 0.9300, 80.0, -0.19),
        :NaNO3 => Salt(i[:Na], 1, i[:NO3], 1, 0.7379, 304.0, -0.39),
        :NH42SO4 => Salt(i[:NH4], 2, i[:SO4], 1, 0.7997, 80.0, -0.25),
        :NH4NO3 => Salt(i[:NH4], 1, i[:NO3], 1, 0.6183, 852.0, -1.15),
        :NH4Cl => Salt(i[:NH4], 1, i[:Cl], 1, 0.7710, 239.0, 0.82),
        :NH4HSO4 => NH4HSO4aqs(i[:NH4], 1, i[:HSO4], 1, 0.4000, 384.0, NaN),
        :NaHSO4 => NaHSO4aqs(i[:Na], 1, i[:HSO4], 1, 0.5200, -45.0, NaN),
        :NH43HSO42 => NH43HSO42aqs(i[:NH4], 3, i[:HSO4], 2, 0.6900, 186.0, Inf),  
        :H2SO4 => Salt(i[:H], 2, i[:SO4], 1, 0.000, missing, -0.1),
        :HHSO4 => Salt(i[:H], 1, i[:HSO4], 1, 0.000, missing, 8.00),
        :HNO3 => Salt(i[:H], 1, i[:NO3], 1, NaN, missing, 2.60), # There is no aqueous to solid conversion for HNO3.
        :HCl => Salt(i[:H], 1, i[:Cl], 1, NaN, missing, 6.00), # There is no aqueous to solid conversion for HCl.
    )
end

"""
Add additional salts as necessary to all the calculation of the activity 
of the SpecialSalts, where s is a dictionary of *all* salts.
"""
function augmented_salts(active_salts, s)
    extra_salts = Dict( # Additional species needed to calculate activity of the SpecialSalts
        s[:KHSO4] => [s[:HHSO4], s[:KCl], s[:HCl]],
        s[:NH4HSO4] => [s[:HHSO4], s[:NH4Cl], s[:HCl]],
        s[:NaHSO4] => [s[:HHSO4], s[:NaCl], s[:HCl]],
        s[:NH43HSO42] => [s[:NH42SO4], s[:NH4HSO4]],
    )
    for es in keys(extra_salts)
        if es in active_salts
            active_salts = vcat(active_salts, extra_salts[es])
        end
    end
    unique(active_salts)
end

"""
Find all salts that have the same cation as the given salt.
"""
function same_cation(s::SaltLike, active_salts, all_salt_dict)
    aug_salts = augmented_salts(active_salts, all_salt_dict)
    aug_salts[[s.cation == ss.cation for ss in aug_salts]]
end

"""
Find all salts that have the same anion as the given salt.
"""
function same_anion(s::SaltLike, active_salts, all_salt_dict)
    aug_salts = augmented_salts(active_salts, all_salt_dict)
    aug_salts[[s.anion == ss.anion for ss in aug_salts]]
end

### Calculate aqueous activity coefficients.

# NOTE: The paper (between equations 6 and 7) says that the units of Aᵧ are kg^0.5 mol^−0.5, but the equation
# doesn't work unless those units are inverted.
@constants Aᵧ = 0.511 [unit = u"mol^0.5/kg_water^0.5", description = "Debye-Hückel constant at 298.15 K"]
@constants I_one = 1 [unit = u"mol/kg_water", description = "An ionic strength of 1"]

# Equation 6
logγ₁₂T⁰(s::Salt, active_salts, sss, I, W) = -Aᵧ * (abs(s.cation.z) * abs(s.anion.z) * √I) / (√I_one + √I) +
                                        (abs(s.cation.z) * abs(s.anion.z)) / (abs(s.cation.z) + abs(s.anion.z)) *
                                        (F₁(s, active_salts, sss, I, W) / abs(s.cation.z) + F₂(s, active_salts, sss, I, W) / abs(s.anion.z))

# Equation 7
F₁(s::Salt, active_salts, sss, I, W) = sum([
    Y(ss, I, W) * logγ⁰₁₂(ss, I) * √I_one +
    Aᵧ * √I / (√I_one + √I) * abs(ss.cation.z) * abs(ss.anion.z) * Y(ss, I, W)
    for ss ∈ same_cation(s, active_salts, sss)
])

# Equation 8
F₂(s::Salt, active_salts, sss, I, W) = sum([
    X(ss, I, W) * logγ⁰₁₂(ss, I) * √I_one +
    Aᵧ * √I / (√I_one + √I) * abs(ss.cation.z) * abs(ss.anion.z) * X(ss, I, W)
    for ss ∈ same_anion(s, active_salts, sss)
])

# Supplemental equations after 7 and 8
Y(s::SaltLike, I, W) = ((abs(s.cation.z) + abs(s.anion.z)) / 2)^2 * s.anion.m / W / I
X(s::SaltLike, I, W) = ((abs(s.cation.z) + abs(s.anion.z)) / 2)^2 * s.cation.m / W / I

# Equation 9
logγ⁰₁₂(s::Salt, I) = abs(s.cation.z) * abs(s.anion.z) * log(Γ⁰(s.q, I) / I_one)
# Equation 10
Γ⁰(q, I) = (I_one + B(q) * ((I_one + 0.1I) / I_one)^q * I_one - I_one * B(q)) * Γstar(q, I)
# Equation 11
B(q) = 0.75 - 0.065q
# Equation 12
Γstar(q, I) = exp(-0.5107√I / (√I_one + C(q, I) * √I))
# Equation 13
C(q, I) = 1 + 0.055q * exp(-0.023I^3 / I_one^3)

@constants T₀₂ = 273.15 [unit = u"K", description = "Standard temperature 2"]
@constants c_1 = 0.005 [unit = u"K^-1", description = "Constant for Fountoukis and Nenes (2007) Eq. 14"]

# Equation 14
A(I) = -((0.41√I / (√I_one + √I)) + 0.039(I / I_one)^0.92)
logγ₁₂(s::Salt, active_salts, sss, I, W) = (1.125 - c_1 * (T - T₀₂)) * logγ₁₂T⁰(s, active_salts, sss, I, W) /
                                      √I_one - (0.125 - c_1 * (T - T₀₂)) * A(I)

"""
Calculate the activity coefficient of a salt based on Section 2.2
in Fountoukis and Nenes (2007).
"""
activity(s::Salt, active_salts, sss, I, W) = (s.cation.m / W)^s.ν_cation * (s.anion.m / W)^s.ν_anion *
                                        exp(logγ₁₂(s, active_salts, sss, I, W))^(s.ν_cation + s.ν_anion)

### Special cases

abstract type SpecialSalt <: SaltLike end

"""
The activity of a SpecialSalt is the same as for a salt except that it has 
a special activity coefficient as defined in the footnotes to Table 4.
"""
activity(s::SpecialSalt, active_salts, sss, I, W) = (s.cation.m / W)^s.ν_cation * (s.anion.m / W)^s.ν_anion *
                                               γ₁₂(s, active_salts, sss, I, W)^(s.ν_cation + s.ν_anion)

logγ⁰₁₂(s::SpecialSalt, I) = abs(s.cation.z) * abs(s.anion.z) * log(Γ⁰(0, I) / I_one) # Assuming q = 0 for SpecialSalts

# Generate special salt structs for multiple dispatch.
specialsaltnames = [:CaSO4, :KHSO4, :NH4HSO4, :NaHSO4, :NH43HSO42]
for (i, name) ∈ enumerate(specialsaltnames)
    structname = Symbol(name, "aqs")
    eval(quote
        """
        From the footnotes to Table 4, CaSO4 has a special activity coefficient.
        """
        struct $structname <: SpecialSalt
            cation::Ion
            "Number of cations per molecule"
            ν_cation::Number
            anion::Ion
            "Number of anions per molecule"
            ν_anion::Number
            "Deliquescence relative humidity at 298.15K"
            drh
            "Enthalpy term (-18/1000R L_s m_s)"
            l_term
            "Kusik-Meissner Binary activity parameter"
            q::Number
        end
    end)
end

"""
From Table 4 footnote a, CaSO4 has an activity coefficient of zero.
"""
γ₁₂(s::CaSO4aqs, active_salts, sss, I, W) = 1.0e-20

"""
From Table 4 footnote b, KHSO4 has a unique activity coefficient
"""
γ₁₂(s::KHSO4aqs, active_salts, sss, I, W) = (exp(logγ₁₂(sss[:HHSO4], active_salts, sss, I, W)) *
                                        exp(logγ₁₂(sss[:KCl], active_salts, sss, I, W)) /
                                        exp(logγ₁₂(sss[:HCl], active_salts, sss, I, W)))^(1 / 2)

"""
From Table 4 footnote c, NH4HSO4 has a unique activity coefficient
"""
γ₁₂(s::NH4HSO4aqs, active_salts, sss, I, W) = (exp(logγ₁₂(sss[:HHSO4], active_salts, sss, I, W)) *
                                          exp(logγ₁₂(sss[:NH4Cl], active_salts, sss, I, W)) /
                                          exp(logγ₁₂(sss[:HCl], active_salts, sss, I, W)))^(1 / 2)

"""
From Table 4 footnote d, NaHSO4 has a unique activity coefficient
"""
γ₁₂(s::NaHSO4aqs, active_salts, sss, I, W) = (exp(logγ₁₂(sss[:HHSO4], active_salts, sss, I, W)) *
                                         exp(logγ₁₂(sss[:NaCl], active_salts, sss, I, W)) /
                                         exp(logγ₁₂(sss[:HCl], active_salts, sss, I, W)))^(1 / 2)

"""
From Table 4 footnote e, NH43HSO42 has a unique activity coefficient
"""
γ₁₂(s::NH43HSO42aqs, active_salts, sss, I, W) = (exp(logγ₁₂(sss[:NH42SO4], active_salts, sss, I, W))^3 *
                                            γ₁₂(sss[:NH4HSO4], active_salts, sss, I, W))^(1 / 5)

"""
Create a system of equations for the ionic strength of the multicomponent solution, as specified by
Fountoukis and Nenes (2007), between equations 8 and 9: ``I = \\frac{1}{2} \\sum_i m_i z_i^2``
"""
function IonicStrength(t, all_Ions, W)
    @variables I(t) = 1.0e-4 [unit = u"mol/kg_water", description = "Ionic strength"]
    I = ParentScope(I)
    @constants I_one = 1 [unit = u"mol/kg_water", description = "An ionic strength of 1"]
    @named ionic_strength = ODESystem([
            I ~ 1 / 2 * sum([ion.m / W * ion.z^2 for ion in all_Ions])
        ], t, [I], [])
end

"""
Create an equation system for the activities of all the species listed in all_species,
where `f` is a function that takes a species as its argument and returns its activity.
"""
function Activities(t, all_species, f)
    eqs = Equation[]
    vars = []
    for s ∈ all_species
        rhs = f(s)
        unit = ModelingToolkit.get_unit(rhs)
        x = Symbol("a_", nameof(s))
        a = only(@variables $x(t) = 1 [unit = unit, description = "Activity of $(nameof(s))"])
        push!(vars, a)
        push!(eqs, a ~ rhs)
    end
    ODESystem(eqs, t, vars, []; name=:activities)
end