"""
An ion with the given concentration `m` and valence (charge) `z`.
"""
struct Ion <: Species
    m::Num # Concentration in mol/(m^3 air)
    m_water # Concentration in mol/(kg water)
    z::Int
end

#==
NOTE: Fountoukis and Nenes (2007) don't give any details on how to calculate the 
activity coefficient of an ion, so we assume that it is 1.0.
==#
γ(i::Ion) = 1.0 / W
terms(i::Ion) = [i.m], [1]

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
        @species $varname($t) = 1e-20 [unit = u"mol/m_air^3", description = $s]
        $varname = ParentScope($varname)
        push!(all_ions, $varname)
        $n = Ion($varname, $varname/W, $(ion_valence[i]))
        push!(all_Ions, $n)
    end)
end

"""
An aqueous salt comprised of a cation, an anion, and an activity parameter (q).
q values are given in Table 4 of Fountoukis and Nenes (2007).
"""
struct Salt <: Species
    cation::Ion
    "Number of cations per molecule"
    ν_cation::Number
    anion::Ion
    "Number of anions per molecule"
    ν_anion::Number
    "Binary activity parameter"
    q::Number

    function Salt(cation::Ion, ν_cation::Number, anion::Ion, ν_anion::Number, q::Number)
        if cation.z * ν_cation + anion.z * ν_anion ≠ 0
            if q != Inf # Special case for NH43HSO42, which doesn't balance.
                throw(ArgumentError("The charge of the cation and anion must sum to zero."))
            end
        end
        new(cation, ν_cation, anion, ν_anion, q)
    end
end

function Base.nameof(s::Salt)
    c = replace(string(Symbolics.tosymbol(s.cation.m, escape=false)), "_aq" => "")
    a = replace(string(Symbolics.tosymbol(s.anion.m, escape=false)), "_aq" => "")
    "$(c)$(s.ν_cation > 1 ? s.ν_cation : "")$(a)$(s.ν_anion > 1 ? s.ν_anion : "")"
end

# Salts from Table 4.
CaNO32_aqs = Salt(Ca_ion, 1, NO3_ion, 2, 0.93)
CaCl2_aqs = Salt(Ca_ion, 1, Cl_ion, 2, 2.4)
# CaSO4 and KHSO4 are below as SpecialSalts.
K2SO4_aqs = Salt(K_ion, 2, SO4_ion, 1, -0.25)
KNO3_aqs = Salt(K_ion, 1, NO3_ion, 1, -2.33)
KCl_aqs = Salt(K_ion, 1, Cl_ion, 1, 0.92)
MgSO4_aqs = Salt(Mg_ion, 1, SO4_ion, 1, 0.15)
MgNO32_aqs = Salt(Mg_ion, 1, NO3_ion, 2, 2.32)
MgCl2_aqs = Salt(Mg_ion, 1, Cl_ion, 2, 2.90)
NaCl_aqs = Salt(Na_ion, 1, Cl_ion, 1, 2.23)
Na2SO4_aqs = Salt(Na_ion, 2, SO4_ion, 1, -0.19)
NaNO3_aqs = Salt(Na_ion, 1, NO3_ion, 1, -0.39)
NH42SO4_aqs = Salt(NH4_ion, 2, SO4_ion, 1, -0.25)
NH4NO3_aqs = Salt(NH4_ion, 1, NO3_ion, 1, -1.15)
NH4Cl_aqs = Salt(NH4_ion, 1, Cl_ion, 1, 0.82)
# NH4HSO4, NaHSO4, and NH43HSO42 are below as SpecialSalts.
H2SO4_aqs = Salt(H_ion, 2, SO4_ion, 1, -0.1)
HHSO4_aqs = Salt(H_ion, 1, HSO4_ion, 1, 8.00)
HNO3_aqs = Salt(H_ion, 1, NO3_ion, 1, 2.60)
HCl_aqs = Salt(H_ion, 1, Cl_ion, 1, 6.00)
all_salts = [CaNO32_aqs, CaCl2_aqs, K2SO4_aqs, KNO3_aqs, KCl_aqs, MgSO4_aqs,
    MgNO32_aqs, MgCl2_aqs, NaCl_aqs, Na2SO4_aqs, NaNO3_aqs, NH42SO4_aqs, NH4NO3_aqs,
    NH4Cl_aqs, H2SO4_aqs, HHSO4_aqs, HNO3_aqs, HCl_aqs]


"""
Find all salts that have the same cation as the given salt.
"""
same_cation(s::Salt) = all_salts[[s.cation == ss.cation for ss in all_salts]]
"""
Find all salts that have the same anion as the given salt.
"""
same_anion(s::Salt) = all_salts[[s.anion == ss.anion for ss in all_salts]]

### Calculate aqueous activity coefficients.

# NOTE: The paper (between equations 6 and 7) says that the units of Aᵧ are kg^0.5 mol^−0.5, but the equation
# doesn't work unless those units are inverted.
@constants Aᵧ = 0.511 [unit = u"mol^0.5/kg_water^0.5", description = "Debye-Hückel constant at 298.15 K"]
@constants I_one = 1 [unit = u"mol/kg_water", description = "An ionic strength of 1"]

# Equation 6
logγ₁₂T⁰(s::Salt) = -Aᵧ * (abs(s.cation.z) * abs(s.anion.z) * √I) / (√I_one + √I) +
                    (abs(s.cation.z) * abs(s.anion.z)) / (abs(s.cation.z) + abs(s.anion.z)) *
                    (F₁(s) / abs(s.cation.z) + F₂(s) / abs(s.anion.z))

# Equation 7
F₁(s::Salt) = sum([
    Y(ss) * logγ⁰₁₂(ss) * √I_one +
    Aᵧ * √I / (√I_one + √I) * abs(ss.cation.z) * abs(ss.anion.z) * Y(ss)
    for ss ∈ same_cation(s)
])


# Equation 8
F₂(s::Salt) = sum([
    X(ss) * logγ⁰₁₂(ss) * √I_one +
    Aᵧ * √I / (√I_one + √I) * abs(ss.cation.z) * abs(ss.anion.z) * X(ss)
    for ss ∈ same_anion(s)
])

# Supplemental equations after 7 and 8
Y(s::Salt) = ((abs(s.cation.z) + abs(s.anion.z)) / 2)^2 * s.anion.m_water / I
X(s::Salt) = ((abs(s.cation.z) + abs(s.anion.z)) / 2)^2 * s.cation.m_water / I

# Equation 9
logγ⁰₁₂(s::Salt) = abs(s.cation.z) * abs(s.anion.z) * log(Γ⁰(s.q) / I_one)
# Equation 10
Γ⁰(q) = (I_one + B(q) * ((I_one + 0.1I) / I_one)^q * I_one - I_one * B(q)) * Γstar(q)
# Equation 11
B(q) = 0.75 - 0.065q
# Equation 12
Γstar(q) = exp(-0.5107√I / (√I_one + C(q) * √I))
# Equation 13
C(q) = 1 + 0.055q * exp(-0.023I^3 / I_one^3)

@constants T₀₂ = 273.15 [unit = u"K", description = "Standard temperature 2"]
@constants c_1 = 0.005 [unit = u"K^-1"] 

# Equation 14
A = -((0.41√I / (√I_one + √I)) + 0.039(I / I_one)^0.92)
logγ₁₂(s::Salt) = (1.125 - c_1 * (T - T₀₂)) * logγ₁₂T⁰(s) / √I_one - (0.125 - c_1 * (T - T₀₂)) * A
"""
Calculate the activity coefficient of a salt based on Section 2.2 
in Fountoukis and Nenes (2007).
"""
γ(s::Salt) = exp(logγ₁₂(s))^(s.ν_cation + s.ν_anion) / W^(s.ν_cation + s.ν_anion)

terms(s::Salt) = [s.cation.m, s.anion.m], [s.ν_cation, s.ν_anion]

### Special cases

abstract type SpecialSalt <: Species end

"""
The activity of a SpecialSalt is the same as for a salt except that it has 
a special activity coefficient as defined in the footnotes to Table 4.
"""
#activity(s::SpecialSalt) = s.cation.m_water^s.ν_cation * s.anion.m_water^s.ν_anion * γ(s)

γ(s::SpecialSalt) = γ₁₂(s)^(s.ν_cation + s.ν_anion) / W^(s.ν_cation + s.ν_anion)

terms(s::SpecialSalt) = [s.cation.m, s.anion.m], [s.ν_cation, s.ν_anion]

specialsaltnames = [:CaSO4, :KHSO4, :NH4HSO4, :NaHSO4, :NH43HSO42]
specialsalts = [Salt(Ca_ion, 1, SO4_ion, 1, NaN), Salt(K_ion, 1, HSO4_ion, 1, NaN),
    Salt(NH4_ion, 1, HSO4_ion, 1, NaN), Salt(Na_ion, 1, HSO4_ion, 1, NaN), Salt(NH4_ion, 3, HSO4_ion, 2, Inf)]

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
        end
        $varname = $structname($(s.cation), $(s.ν_cation), $(s.anion), $(s.ν_anion))
    end)
end

"""
From Table 4 footnote a, CaSO4 has an activity coefficient of zero.
"""
γ₁₂(s::CaSO4aqs) = 1.0e-20

"""
From Table 4 footnote b, KHSO4 has a unique activity coefficient
"""
γ₁₂(s::KHSO4aqs) = (exp(logγ₁₂(HHSO4_aqs)) * exp(logγ₁₂(KCl_aqs)) / exp(logγ₁₂(HCl_aqs)))^(1 / 2)

"""
From Table 4 footnote c, NH4HSO4 has a unique activity coefficient
"""
γ₁₂(s::NH4HSO4aqs) = (exp(logγ₁₂(HHSO4_aqs)) * exp(logγ₁₂(NH4Cl_aqs)) / exp(logγ₁₂(HCl_aqs)))^(1 / 2)

"""
From Table 4 footnote d, NaHSO4 has a unique activity coefficient
"""
γ₁₂(s::NaHSO4aqs) = (exp(logγ₁₂(HHSO4_aqs)) * exp(logγ₁₂(NaCl_aqs)) / exp(logγ₁₂(HCl_aqs)))^(1 / 2)

"""
From Table 4 footnote e, NH43HSO42 has a unique activity coefficient
"""
γ₁₂(s::NH43HSO42aqs) = (exp(logγ₁₂(NH42SO4_aqs))^3 * γ₁₂(NH4HSO4_aqs))^(1 / 5)


# Tests 
@test length(all_ions) == 14
@test length(all_salts) == 18

@test same_cation(KCl_aqs) == [K2SO4_aqs, KNO3_aqs, KCl_aqs]
@test same_anion(KCl_aqs) == [CaCl2_aqs, KCl_aqs, MgCl2_aqs, NaCl_aqs, NH4Cl_aqs, HCl_aqs]

@test ModelingToolkit.get_unit(logγ₁₂T⁰(NaCl_aqs)) == u"mol^0.5/kg_water^0.5"
@test ModelingToolkit.get_unit(F₁(NaCl_aqs)) == u"mol^0.5/kg_water^0.5"
@test ModelingToolkit.get_unit(F₂(NaCl_aqs)) == u"mol^0.5/kg_water^0.5"
@test ModelingToolkit.get_unit(X(NaCl_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(Y(NaCl_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(Γ⁰(NaCl_aqs.q)) == u"mol/kg_water"
@test ModelingToolkit.get_unit(Γstar(1)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(C(1)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(logγ⁰₁₂(NaCl_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(logγ⁰₁₂(NaCl_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}
@test ModelingToolkit.get_unit(logγ₁₂(NaCl_aqs)) isa Unitful.FreeUnits{(),NoDims,nothing}

# Activity coefficients should be ≤ 1 for ions in aqueous solution
@test substitute(logγ⁰₁₂(KCl_aqs), [I => 2.5, Aᵧ => 0.511, Cl_aq => 1, I_one => 1, H_aq => 2,
    K_aq => 0.5, Mg_aq => 1.2, NH4_aq => 2.5, NO3_aq => 3.1, Na_aq => 0.2, Ca_aq => 1.2,
    SO4_aq => 2.0]) ≈ -0.16013845145909214

@test substitute(logγ₁₂T⁰(KCl_aqs), [I => 2.5, Aᵧ => 0.511, Cl_aq => 1, I_one => 1, H_aq => 2,
    K_aq => 0.5, Mg_aq => 1.2, NH4_aq => 2.5, NO3_aq => 3.1, Na_aq => 0.2, Ca_aq => 1.2, W => 1,
    SO4_aq => 2.0]) ≈ 0.6031282826436271

# Activity coefficients should decrease with increasing temperature.
@test substitute(logγ₁₂(KCl_aqs), [I => 2.5, Aᵧ => 0.511, Cl_aq => 1, I_one => 1, H_aq => 2,
    K_aq => 0.5, Mg_aq => 1.2, NH4_aq => 2.5, NO3_aq => 3.1, Na_aq => 0.2, Ca_aq => 1.2, W => 1,
    SO4_aq => 2.0, T => 310, T₀₂ => 273.15, c_1 => 0.005]) ≈ 0.5471434157637021

@test substitute(logγ₁₂(KCl_aqs), [I => 2.5, Aᵧ => 0.511, Cl_aq => 1, I_one => 1, H_aq => 2,
    K_aq => 0.5, Mg_aq => 1.2, NH4_aq => 2.5, NO3_aq => 3.1, Na_aq => 0.2, Ca_aq => 1.2, W => 1,
    SO4_aq => 2.0, T => 280, T₀₂ => 273.15, c_1 => 0.005]) ≈ 0.6888772559660434

# Activity coefficients should decrease with increasing ionic strength.
@test substitute(logγ₁₂(KCl_aqs), [I => 5, Aᵧ => 0.511, Cl_aq => 1, I_one => 1, H_aq => 2,
    K_aq => 0.5, Mg_aq => 1.2, NH4_aq => 2.5, NO3_aq => 3.1, Na_aq => 0.2, Ca_aq => 1.2, W => 1,
    SO4_aq => 2.0, T => 298.15, T₀₂ => 273.15, c_1 => 0.005]) ≈ 0.5664133683846955

# Test activity coefficients for all salts.
want_vals = [0.4477168981554961, 1.3684551018017284, -1.125788329292307, -0.12414028434997845, 
    0.5664133683846955, -0.022779782656542125, 0.9450814283879742, 1.865819632034206, 0.831071070150315, 
    -0.772911393604814, 0.140517417415641, -0.9501106424586302, 0.007617980775278976, 0.6981716335099528, 
    -0.008228428441089064, 0.5751777380058343, 0.7140296412884348, 1.4045832940231087]
for (i, salt) in enumerate(all_salts)
    v = substitute(logγ₁₂(salt), [I => 5, Aᵧ => 0.511, Cl_aq => 1, I_one => 1, H_aq => 2, W => 1,
        K_aq => 0.5, Mg_aq => 1.2, NH4_aq => 2.5, NO3_aq => 3.1, Na_aq => 0.2, Ca_aq => 1.2,
        SO4_aq => 2.0, HSO4_aq => 1.0, T => 298.15, T₀₂ => 273.15, c_1 => 0.005])
    @test v ≈ want_vals[i]
end

# Units in last column in Table 2.
@test ModelingToolkit.get_unit(activity(NaCl_aqs)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(activity(CaNO32_aqs)) == u"mol^3/kg_water^3"

@test substitute(activity(NaCl_aqs), [I => 5, Aᵧ => 0.511, Cl_aq => 1, I_one => 1, H_aq => 2,
    K_aq => 0.5, Mg_aq => 1.2, NH4_aq => 2.5, NO3_aq => 3.1, Na_aq => 0.2, Ca_aq => 1.2, W => 1,
    SO4_aq => 2.0, HSO4_aq => 1.0, T => 298.15, T₀₂ => 273.15, c_1 => 0.005]) ≈ 1.0541178203304393

# Special cases

# Units in last column in Table 2.
@test ModelingToolkit.get_unit(activity(CaSO4_aqs)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(activity(KHSO4_aqs)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(activity(NH4HSO4_aqs)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(activity(NaHSO4_aqs)) == u"mol^2/kg_water^2"
@test ModelingToolkit.get_unit(activity(NH43HSO42_aqs)) == u"mol^5/kg_water^5"

want_γ = [1.0e-20, 0.876782700005408, 0.9364895073669987, 1.000833103904957, 0.5581152605946901]
want_activity = [2.4e-40, 0.3843739515143867, 2.192531493521209, 0.20033338037440607, 0.8461347162329317]
for (i, s) ∈ enumerate([CaSO4_aqs, KHSO4_aqs, NH4HSO4_aqs, NaHSO4_aqs, NH43HSO42_aqs])
    @test substitute(γ₁₂(s), [I => 5, Aᵧ => 0.511, Cl_aq => 1, I_one => 1, H_aq => 2, W => 1,
        K_aq => 0.5, Mg_aq => 1.2, NH4_aq => 2.5, NO3_aq => 3.1, Na_aq => 0.2, Ca_aq => 1.2,
        SO4_aq => 2.0, HSO4_aq => 1.0, T => 298.15, T₀₂ => 273.15, c_1 => 0.005]) ≈ want_γ[i]

    @test substitute(activity(s), [I => 5, Aᵧ => 0.511, Cl_aq => 1, I_one => 1, H_aq => 2, W => 1,
        K_aq => 0.5, Mg_aq => 1.2, NH4_aq => 2.5, NO3_aq => 3.1, Na_aq => 0.2, Ca_aq => 1.2,
        SO4_aq => 2.0, HSO4_aq => 1.0, T => 298.15, T₀₂ => 273.15, c_1 => 0.005]) ≈ want_activity[i]
end

@test nameof(CaCl2_aqs) == "CaCl2"