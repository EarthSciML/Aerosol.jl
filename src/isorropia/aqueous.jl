"""
An ion with the given concentration `m` and valence (charge) `z`.
"""
struct Ion <: Species
    m::Num # Concentration in mol/(kg water)
    z::Int
end

function Base.nameof(i::Ion)
    string(Symbolics.tosymbol(i.m, escape=false))
end

vars(i::Ion) = [i.m]

#==
The activity coefficient of an ion is assumed to be one (Fountoukis and Nenes (2007), Section 3.3).
==#
activity(i::Ion) = i.m

# Generate the aqueous ions.
# Each ion has an associated MTK variable named 
# <name>_aq, where <name> is the name of the compound, and
# a Ion struct named <name>_ion.
ion_names = [:NH4, :Na, :H, :Cl, :NO3, :SO4, :HNO3, :NH3, :HCl, :HSO4, :Ca, :K, :Mg, :OH]
ion_valence = [1, 1, 1, -1, -1, -2, 0, 0, 0, -1, 2, 1, 2, -1]
ion_charge = [1, 1, 1, -1, -1, -2, 0,]
all_ions = []
all_Ions = []
for i ∈ eachindex(ion_names)
    n = Symbol(ion_names[i], "_ion")
    varname = Symbol(ion_names[i], "_aq")
    s = "Aqueous $(ion_names[i])"
    eval(quote
        @species $varname = 1e-11 [unit = u"mol/kg_water", description = $s]
        $varname = ParentScope($varname)
        push!(all_ions, $varname)
        $n = Ion($varname, $(ion_valence[i]))
        push!(all_Ions, $n)
    end)
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

# Salts from Table 4.
CaNO32_aqs = Salt(Ca_ion, 1, NO3_ion, 2, 0.4906, 509.4, 0.93)
CaCl2_aqs = Salt(Ca_ion, 1, Cl_ion, 2, 0.2830, 551.1, 2.4)
# CaSO4 and KHSO4 are below as SpecialSalts.
K2SO4_aqs = Salt(K_ion, 2, SO4_ion, 1, 0.9751, 35.6, -0.25)
KNO3_aqs = Salt(K_ion, 1, NO3_ion, 1, 0.9248, missing, -2.33)
KCl_aqs = Salt(K_ion, 1, Cl_ion, 1, 0.8426, 158.9, 0.92)
MgSO4_aqs = Salt(Mg_ion, 1, SO4_ion, 1, 0.8613, -714.5, 0.15)
MgNO32_aqs = Salt(Mg_ion, 1, NO3_ion, 2, 0.5400, 230.2, 2.32)
MgCl2_aqs = Salt(Mg_ion, 1, Cl_ion, 2, 0.3284, 42.23, 2.90)
NaCl_aqs = Salt(Na_ion, 1, Cl_ion, 1, 0.7528, 25.0, 2.23)
Na2SO4_aqs = Salt(Na_ion, 2, SO4_ion, 1, 0.9300, 80.0, -0.19)
NaNO3_aqs = Salt(Na_ion, 1, NO3_ion, 1, 0.7379, 304.0, -0.39)
NH42SO4_aqs = Salt(NH4_ion, 2, SO4_ion, 1, 0.7997, 80.0, -0.25)
NH4NO3_aqs = Salt(NH4_ion, 1, NO3_ion, 1, 0.6183, 852.0, -1.15)
NH4Cl_aqs = Salt(NH4_ion, 1, Cl_ion, 1, 0.7710, 239.0, 0.82)
# NH4HSO4, NaHSO4, and NH43HSO42 are below as SpecialSalts.
H2SO4_aqs = Salt(H_ion, 2, SO4_ion, 1, 0.000, missing, -0.1)
HHSO4_aqs = Salt(H_ion, 1, HSO4_ion, 1, 0.000, missing, 8.00)
HNO3_aqs = Salt(H_ion, 1, NO3_ion, 1, NaN, missing, 2.60) # There is no aqueous to solid conversion for HNO3.
HCl_aqs = Salt(H_ion, 1, Cl_ion, 1, NaN, missing, 6.00) # There is no aqueous to solid conversion for HCl.
all_salts = SaltLike[CaNO32_aqs, CaCl2_aqs, K2SO4_aqs, KNO3_aqs, KCl_aqs, MgSO4_aqs,
    MgNO32_aqs, MgCl2_aqs, NaCl_aqs, Na2SO4_aqs, NaNO3_aqs, NH42SO4_aqs, NH4NO3_aqs,
    NH4Cl_aqs, H2SO4_aqs, HHSO4_aqs, HNO3_aqs, HCl_aqs]

"""
Add additional salts as necessary to all the calculation of the activity 
of the SpecialSalts.
"""
function augmented_salts(all_salts)
    extra_salts = Dict( # Additional species needed to calculate activity of the SpecialSalts
        KHSO4_aqs => [HHSO4_aqs, KCl_aqs, HCl_aqs],
        NH4HSO4_aqs => [HHSO4_aqs, NH4Cl_aqs, HCl_aqs],
        NaHSO4_aqs => [HHSO4_aqs, NaCl_aqs, HCl_aqs],
        NH43HSO42_aqs => [NH42SO4_aqs, NH4HSO4_aqs],
    )
    for es in keys(extra_salts)
        if es in all_salts
            all_salts = vcat(all_salts, extra_salts[es])
        end
    end
    unique(all_salts)
end

"""
Find all salts that have the same cation as the given salt.
"""
function same_cation(s::SaltLike, all_salts)
    aug_salts = augmented_salts(all_salts)
    aug_salts[[s.cation == ss.cation for ss in aug_salts]]
end

"""
Find all salts that have the same anion as the given salt.
"""
function same_anion(s::SaltLike, all_salts)
    aug_salts = augmented_salts(all_salts)
    aug_salts[[s.anion == ss.anion for ss in aug_salts]]
end

### Calculate aqueous activity coefficients.

# NOTE: The paper (between equations 6 and 7) says that the units of Aᵧ are kg^0.5 mol^−0.5, but the equation
# doesn't work unless those units are inverted.
@constants Aᵧ = 0.511 [unit = u"mol^0.5/kg_water^0.5", description = "Debye-Hückel constant at 298.15 K"]
@constants I_one = 1 [unit = u"mol/kg_water", description = "An ionic strength of 1"]

# Equation 6
logγ₁₂T⁰(s::Salt, all_salts, I) = -Aᵧ * (abs(s.cation.z) * abs(s.anion.z) * √I) / (√I_one + √I) +
                                  (abs(s.cation.z) * abs(s.anion.z)) / (abs(s.cation.z) + abs(s.anion.z)) *
                                  (F₁(s, all_salts, I) / abs(s.cation.z) + F₂(s, all_salts, I) / abs(s.anion.z))

# Equation 7
F₁(s::Salt, all_salts, I) = sum([
    Y(ss, I) * logγ⁰₁₂(ss, I) * √I_one +
    Aᵧ * √I / (√I_one + √I) * abs(ss.cation.z) * abs(ss.anion.z) * Y(ss, I)
    for ss ∈ same_cation(s, all_salts)
])


# Equation 8
F₂(s::Salt, all_salts, I) = sum([
    X(ss, I) * logγ⁰₁₂(ss, I) * √I_one +
    Aᵧ * √I / (√I_one + √I) * abs(ss.cation.z) * abs(ss.anion.z) * X(ss, I)
    for ss ∈ same_anion(s, all_salts)
])

# Supplemental equations after 7 and 8
Y(s::SaltLike, I) = ((abs(s.cation.z) + abs(s.anion.z)) / 2)^2 * s.anion.m / I
X(s::SaltLike, I) = ((abs(s.cation.z) + abs(s.anion.z)) / 2)^2 * s.cation.m / I

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
logγ₁₂(s::Salt, all_salts, I) = (1.125 - c_1 * (T - T₀₂)) * logγ₁₂T⁰(s, all_salts, I) / √I_one - (0.125 - c_1 * (T - T₀₂)) * A(I)

"""
Calculate the activity coefficient of a salt based on Section 2.2
in Fountoukis and Nenes (2007).
"""
activity(s::Salt, all_salts, I) = s.cation.m^s.ν_cation * s.anion.m^s.ν_anion #*
                                  #exp(logγ₁₂(s, all_salts, I))^(s.ν_cation + s.ν_anion)

### Special cases

abstract type SpecialSalt <: SaltLike end

"""
The activity of a SpecialSalt is the same as for a salt except that it has 
a special activity coefficient as defined in the footnotes to Table 4.
"""
activity(s::SpecialSalt, all_salts, I) = s.cation.m^s.ν_cation * s.anion.m^s.ν_anion #*
                                         #γ₁₂(s, all_salts, I)^(s.ν_cation + s.ν_anion)
logγ⁰₁₂(s::SpecialSalt, I) = abs(s.cation.z) * abs(s.anion.z) * log(Γ⁰(0, I) / I_one) # Assuming q = 0 for SpecialSalts

min_conc(s::SpecialSalt) = min(s.cation.m, s.anion.m)

specialsaltnames = [:CaSO4, :KHSO4, :NH4HSO4, :NaHSO4, :NH43HSO42]
# Data from Fountoukis and Nenes (2007) Table 4
specialsalts = [
    Salt(Ca_ion, 1, SO4_ion, 1, 0.9700, missing, NaN),
    Salt(K_ion, 1, HSO4_ion, 1, 0.8600, missing, NaN),
    Salt(NH4_ion, 1, HSO4_ion, 1, 0.4000, 384.0, NaN),
    Salt(Na_ion, 1, HSO4_ion, 1, 0.5200, -45.0, NaN),
    Salt(NH4_ion, 3, HSO4_ion, 2, 0.6900, 186.0, Inf)]

for (i, name) ∈ enumerate(specialsaltnames)
    s = specialsalts[i]
    varname = Symbol(name, "_aqs")
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
        end
        $varname = $structname($(s.cation), $(s.ν_cation), $(s.anion), $(s.ν_anion), $(s.drh), $(s.l_term))
        push!(all_salts, $varname)
    end)
end

"""
From Table 4 footnote a, CaSO4 has an activity coefficient of zero.
"""
γ₁₂(s::CaSO4aqs, all_salts, I) = 1.0e-20

"""
From Table 4 footnote b, KHSO4 has a unique activity coefficient
"""
γ₁₂(s::KHSO4aqs, all_salts, I) = (exp(logγ₁₂(HHSO4_aqs, all_salts, I)) * exp(logγ₁₂(KCl_aqs, all_salts, I)) /
                                  exp(logγ₁₂(HCl_aqs, all_salts, I)))^(1 / 2)

"""
From Table 4 footnote c, NH4HSO4 has a unique activity coefficient
"""
γ₁₂(s::NH4HSO4aqs, all_salts, I) = (exp(logγ₁₂(HHSO4_aqs, all_salts, I)) * exp(logγ₁₂(NH4Cl_aqs, all_salts, I)) /
                                    exp(logγ₁₂(HCl_aqs, all_salts, I)))^(1 / 2)

"""
From Table 4 footnote d, NaHSO4 has a unique activity coefficient
"""
γ₁₂(s::NaHSO4aqs, all_salts, I) = (exp(logγ₁₂(HHSO4_aqs, all_salts, I)) * exp(logγ₁₂(NaCl_aqs, all_salts, I)) /
                                   exp(logγ₁₂(HCl_aqs, all_salts, I)))^(1 / 2)

"""
From Table 4 footnote e, NH43HSO42 has a unique activity coefficient
"""
γ₁₂(s::NH43HSO42aqs, all_salts, I) = (exp(logγ₁₂(NH42SO4_aqs, all_salts, I))^3 * γ₁₂(NH4HSO4_aqs, all_salts, I))^(1 / 5)

"""
Create a system of equations for the ionic strength of the multicomponent solution, as specified by
Fountoukis and Nenes (2007), between equations 8 and 9: ``I = \\frac{1}{2} \\sum_i m_i z_i^2``
"""
function IonicStrength(all_Ions)
    @variables I = 1.0e-4 [unit = u"mol/kg_water", description = "Ionic strength"]
    I = ParentScope(I)
    @constants I_one = 1 [unit = u"mol/kg_water", description = "An ionic strength of 1"]
    @named ionic_strength = NonlinearSystem([
            # Force I to always be positive to avoid attempts to take the square root of a negative number.
            I ~ max(1.0e-20 * I_one, 1 / 2 * sum([ion.m * ion.z^2 for ion in all_Ions]))
        ], [I], [I_one])
end

"""
Create an equation system for the activities of all the salts listed in all_salts,
using I as the ionic strength.
"""
function SaltActivities(all_salts, I)
    eqs = Equation[]
    vars = []
    for s ∈ all_salts
        rhs = activity(s, all_salts, I)
        unit = ModelingToolkit.get_unit(rhs)
        x = Symbol("a_", nameof(s))
        a = only(@variables $x=1 [unit = unit, description = "Activity of $(nameof(s))"])
        push!(vars, a)
        push!(eqs, a ~ activity(s, all_salts, I))
    end
    NonlinearSystem(eqs, vars, [I_one, T₀₂, c_1, Aᵧ]; name=:activities)
end


"""
Create an equation system for the activities of all the species listed in all_species.
"""
function Activities(all_species)
    eqs = Equation[]
    vars = []
    for s ∈ all_species
        rhs = activity(s)
        unit = ModelingToolkit.get_unit(rhs)
        x = Symbol("a_", nameof(s))
        a = only(@variables $x=1 [unit = unit, description = "Activity of $(nameof(s))"])
        push!(vars, a)
        push!(eqs, a ~ rhs)
    end
    NonlinearSystem(eqs, vars, []; name=:activities)
end