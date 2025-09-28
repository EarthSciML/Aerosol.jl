using ModelingToolkit
using ModelingToolkit: t
using DynamicQuantities

@mtkmodel Ion begin
    @description "An aqueous ion."
    @parameters begin
        z, [description = "Valence (charge) of the ion"]
    end
    @variables begin
        m(t) = 0.0, [description = "Concentration of the ion", unit = u"mol/m^3"]
        activity(t),
        [
            description = "Activity of the ion. The activity coefficient of an ion is assumed to be one (Fountoukis and Nenes (2007), Section 3.3).",
            unit = u"mol/m^3"]
    end
    @equations begin
        activity ~ m
    end
end

@mtkmodel Salt begin
    @description """
An aqueous salt comprised of a cation, an anion, and an activity parameter (q).
q values are given in Table 4 of Fountoukis and Nenes (2007).
"""
    @structural_parameters begin
        cation = Ion(; z)
        anion = Ion(; z)
    end
    @constants begin
        # Salt properties
        ν_c, [description = "Number of cations per molecule"]
        ν_a, [description = "Number of anions per molecule"]
        drh, [description = "Deliquescence relative humidity at 298.15K"]
        l_t, [description = "Enthalpy term (-18/1000R L_s m_s)"]
        q, [description = "Kusik-Meissner Binary activity parameter"]

        # Derived constants
        logγ⁰, [description = "Log of the standard state activity coefficient"]
        Γ⁰, [unit = u"mol/kg"]
        B
        C
        zz, [description = "Product of the absolute values of the cation and anion charges"]

        # Unit conversions
        I_one = 1, [unit = u"mol/kg", description = "An ionic strength of 1"]
    end
    @variables begin
        X(t), [unit = u"kg/m^3"]
        Y(t), [unit = u"kg/m^3"]
        I(t), [description = "Ionic strength", unit = u"mol/kg"]
        Γ⁺(t)
    end
    @equations begin
        zz ~ ParentScope(cation.z) * ParentScope(anion.z)
        # Supplemental equations after equations 7 and 8
        Y ~ ((ParentScope(cation.z) + ParentScope(anion.z)) / 2)^2 * ParentScope(anion.m) / I
        X ~ ((ParentScope(cation.z) + ParentScope(anion.z)) / 2)^2 * ParentScope(cation.m) / I
        # Equation 9
        logγ⁰ ~ zz * log(Γ⁰ / I_one)
        # Equation 10
        Γ⁰ ~ (I_one + B * ((I_one + 0.1I) / I_one)^q * I_one - I_one * B) * Γ⁺ # TODO(CT): Maybe this should be logΓ⁰ ~ (I_one + B * ((I_one + 0.1I) / I_one)^q * I_one - I_one * B) * logΓ⁺ ?
        # Equation 11
        B ~ 0.75 - 0.065q
        # Equation 12
        Γ⁺ ~ exp(-0.5107√I / (√I_one + C * √I))
        # Equation 13
        C ~ 1 + 0.055q * exp(-0.023I^3 / I_one^3)
    end
end

@mtkmodel Aqueous begin
    @description "Aqueous behavior"
    @constants begin
        # NOTE: The paper (between equations 6 and 7) says that the units of Aᵧ are kg^0.5 mol^−0.5, but the equation
        # doesn't work unless those units are inverted, and Bromley (1973) agrees with that.
        Aᵧ = 0.511,
        [
            unit = u"mol^0.5/kg^0.5",
            description = "Debye-Hückel constant at 298.15 K"
        ]
        I_one = 1, [unit = u"mol/kg", description = "An ionic strength of 1"]
    end
    @variables begin
        I(t), [description = "Ionic strength", unit = u"mol/kg"]

        #! format: off
        F_Ca(t), [description = "Activity contribution from Ca-containing salts", unit = u"kg/m^3"]
        #! format: on
        A_γterm(t), [description = "Debye-Hückel term used in Equation 7 and 8"]
    end
    @components begin
        # Cations
        NH4 = Ion(z = 1)
        Na = Ion(z = 1)
        H = Ion(z = 1)
        Ca = Ion(z = 2)
        K = Ion(z = 1)
        Mg = Ion(z = 2)

        # Anions
        Cl = Ion(z = abs(-1))
        NO3 = Ion(z = abs(-1))
        SO4 = Ion(z = abs(-2))
        HNO3 = Ion(z = 0)
        NH3 = Ion(z = 0)
        HCl = Ion(z = 0)
        HSO4 = Ion(z = abs(-1))
        OH = Ion(z = abs(-1))

        # Salts
        #! format: off
        CaNO32 = Salt(cation = Ca, ν_c = 1, anion = NO3, ν_a = 2, drh = 0.4906, l_t = 509.4, q = 0.93, I=I)
        CaCl2 =  Salt(cation = Ca, ν_c = 1, anion = Cl, ν_a = 2, drh = 0.2830, l_t = 551.1, q = 2.4, I=I)
        CaSO4 = Salt(cation = Ca, ν_c = 1, anion = SO4, ν_a = 1, drh = 0.9700, l_t = missing, q = missing, I=I)
        KHSO4 = Salt(cation = K, ν_c = 1, anion = HSO4, ν_a = 1, drh = 0.8600, l_t = missing, q = missing, I=I)
        K2SO4 = Salt(cation = K, ν_c = 2, anion = SO4, ν_a = 1, drh = 0.9751, l_t = 35.6, q = -0.25, I=I)
        KNO3 = Salt(cation = K, ν_c = 1, anion = NO3, ν_a = 1, drh = 0.9248, l_t = missing, q = -2.33, I=I)
        KCl = Salt(cation = K, ν_c = 1, anion = Cl, ν_a = 1, drh = 0.8426, l_t = 158.9, q = 0.92, I=I)
        MgSO4 = Salt(cation = Mg, ν_c = 1, anion = SO4, ν_a = 1, drh = 0.8613, l_t = -714.5, q = 0.15, I=I)
        MgNO32 = Salt(cation = Mg, ν_c = 1, anion = NO3, ν_a = 2, drh = 0.5400, l_t = 230.2, q = 2.32, I=I)
        MgCl2 = Salt(cation = Mg, ν_c = 1, anion = Cl, ν_a = 2, drh = 0.3284, l_t = 42.23, q = 2.90, I=I)
        NaCl = Salt(cation = Na, ν_c = 1, anion = Cl, ν_a = 1, drh = 0.7528, l_t = 25.0, q = 2.23, I=I)
        Na2SO4 = Salt(cation = Na, ν_c = 2, anion = SO4, ν_a = 1, drh = 0.9300, l_t = 80.0, q = -0.19, I=I)
        NaNO3 = Salt(cation = Na, ν_c = 1, anion = NO3, ν_a = 1, drh = 0.7379, l_t = 304.0, q = -0.39, I=I)
        NH42SO4 = Salt(cation = NH4, ν_c = 2, anion = SO4, ν_a = 1, drh = 0.7997, l_t = 80.0, q = -0.25, I=I)
        NH4NO3 = Salt(cation = NH4, ν_c = 1, anion = NO3, ν_a = 1, drh = 0.6183, l_t = 852.0, q = -1.15, I=I)
        NH4Cl = Salt(cation = NH4, ν_c = 1, anion = Cl, ν_a = 1, drh = 0.7710, l_t = 239.0, q = 0.82, I=I)
        NH4HSO4 = Salt(cation = NH4, ν_c = 1, anion = HSO4, ν_a = 1, drh = 0.4000, l_t = 384.0, q = missing, I=I)
        NaHSO4 = Salt(cation = Na, ν_c = 1, anion = HSO4, ν_a = 1, drh = 0.5200, l_t = -45.0, q = missing, I=I)
        NH43HSO42 = Salt(cation = NH4, ν_c = 3, anion = HSO4, ν_a = 2, drh = 0.6900, l_t = 186.0, q = missing, I=I)
        H2SO4 = Salt(cation = H, ν_c = 2, anion = SO4, ν_a = 1, drh = 0.000, l_t = missing, q = -0.1, I=I)
        HHSO4 = Salt(cation = H, ν_c = 1, anion = HSO4, ν_a = 1, drh = 0.000, l_t = missing, q = 8.00, I=I)
        HNO3_aqs = Salt(cation = H, ν_c = 1, anion = NO3, ν_a = 1, drh = NaN, l_t = missing, q = 2.60, I=I) # There is no aqueous to solid conversion for HNO3.
        HCl_aqs = Salt(cation = H, ν_c = 1, anion = Cl, ν_a = 1, drh = NaN, l_t = missing, q = 6.00, I=I) # There is no aqueous to solid conversion for HCl.
        #! format: on
    end
    @equations begin
        A_γterm ~ Aᵧ * √I / (√I_one + √I) / √I_one
        # Equation 7
        F_Ca ~ CaNO32.Y * CaNO32.logγ⁰ + CaCl2.Y * CaCl2.logγ⁰ + CaSO4.Y * CaSO4.logγ⁰ +
               A_γterm * (CaNO32.zz * CaNO32.Y + CaCl2.zz * CaCl2.Y + CaSO4.zz * CaSO4.Y)
    end
end

@named aq = Aqueous()
equations(aq)
aq.OH

structural_simplify(aq)

F_Ca = sum([s[3].Y * s[3].logγ⁰ for s in filter(s -> s[1] == Ca, salts)]) +
       sum([s[3].zz * s[3].logγ⁰ for s in filter(s -> s[1] == Ca, salts)])

salts = [
    (Ca, NO3, CaNO32),
    (Ca, Cl, CaCl2),
    (Ca, SO4, CaSO4),
    (K, HSO4, KHSO4),
    (K, SO4, K2SO4),
    (K, NO3, KNO3),
    (K, Cl, KCl),
    (Mg, SO4, MgSO4),
    (Mg, NO3, MgNO32),
    (Mg, Cl, MgCl2),
    (Na, Cl, NaCl),
    (Na, SO4, Na2SO4),
    (Na, NO3, NaNO3),
    (NH4, SO4, NH42SO4),
    (NH4, NO3, NH4NO3),
    (NH4, Cl, NH4Cl),
    (NH4, HSO4, NH4HSO4),
    (Na, HSO4, NaHSO4),
    (NH4, HSO4, NH43HSO42),
    (H, SO4, H2SO4),
    (H, HSO4, HHSO4),
    (H, NO3, HNO3_aqs),
    (H, Cl, HCl_aqs)
]

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
    drh::Any
    "Enthalpy term (-18/1000R L_s m_s)"
    l_term::Any
    "Kusik-Meissner Binary activity parameter"
    q::Number

    function Salt(
            cation::Ion,
            ν_cation::Number,
            anion::Ion,
            ν_anion::Number,
            drh,
            l_term,
            q::Number
    )
        if cation.z * ν_cation + anion.z * ν_anion ≠ 0
            if q != Inf # Special case for NH43HSO42, which doesn't balance.
                throw(ArgumentError("The charge of the cation and anion must sum to zero."))
            end
        end
        new(cation, ν_cation, anion, ν_anion, drh, l_term, q)
    end
end

# function Base.nameof(s::SaltLike)
#     c = replace(string(Symbolics.tosymbol(s.cation.m, escape = false)), "_aq" => "")
#     a = replace(string(Symbolics.tosymbol(s.anion.m, escape = false)), "_aq" => "")
#     "$(c)$(s.ν_cation > 1 ? s.ν_cation : "")$(a)$(s.ν_anion > 1 ? s.ν_anion : "")_aqs"
# end

# vars(s::SaltLike) = [s.cation.m, s.anion.m]
# terms(s::SaltLike) = [s.cation.m, s.anion.m], [s.ν_cation, s.ν_anion]

"""
Generate Salts from Table 4, where `i` is a dictionary of ions.
"""
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
        :HCl => Salt(i[:H], 1, i[:Cl], 1, NaN, missing, 6.00) # There is no aqueous to solid conversion for HCl.
    )
end

# function get_cation(salt)
#     Dict(
#         CaNO32 => :Ca,
#         CaCl2 => :Ca,
#         CaSO4 => :Ca,
#         KHSO4 => :K,
#         K2SO4 => :K,
#         KNO3 => :K,
#         KCl => :K,
#         MgSO4 => :Mg,
#         MgNO32 => :Mg,
#         MgCl2 => :Mg,
#         NaCl => :Na,
#         Na2SO4 => :Na,
#         NaNO3 => :Na,
#         NH42SO4 => :NH4,
#         NH4NO3 => :NH4,
#         NH4Cl => :NH4,
#         NH4HSO4 => :NH4,
#         NaHSO4 => :Na,
#         NH43HSO42 => :NH4,
#         H2SO4 => :H,
#         HHSO4 => :H,
#         HNO3 => :H,
#         HCl => :H
#     )[salt]
# end

# function get_anion(salt)
#     Dict(
#         CaNO32 => :NO3,
#         CaCl2 => :Cl,
#         CaSO4 => :SO4,
#         KHSO4 => :HSO4,
#         K2SO4 => :SO4,
#         KNO3 => :NO3,
#         KCl => :Cl,
#         MgSO4 => :SO4,
#         MgNO32 => :NO3,
#         MgCl2 => :Cl,
#         NaCl => :Cl,
#         Na2SO4 => :SO4,
#         NaNO3 => :NO3,
#         NH42SO4 => :SO4,
#         NH4NO3 => :NO3,
#         NH4Cl => :Cl,
#         NH4HSO4 => :HSO4,
#         NaHSO4 => :HSO4,
#         NH43HSO42 => :HSO4,
#         H2SO4 => :SO4,
#         HHSO4 => :HSO4,
#         HNO3 => :NO3,
#         HCl => :Cl
#     )[salt]
# end

"""
Add additional salts as necessary to all the calculation of the activity
of the SpecialSalts, where s is a dictionary of *all* salts.
"""
function augmented_salts(active_salts, s)
    extra_salts = Dict( # Additional species needed to calculate activity of the SpecialSalts
        s[:KHSO4] => [s[:HHSO4], s[:KCl], s[:HCl]],
        s[:NH4HSO4] => [s[:HHSO4], s[:NH4Cl], s[:HCl]],
        s[:NaHSO4] => [s[:HHSO4], s[:NaCl], s[:HCl]],
        s[:NH43HSO42] => [s[:NH42SO4], s[:NH4HSO4]]
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

# Equation 6
function logγ₁₂T⁰(s::Salt, active_salts, sss, I, W)
    -Aᵧ * (abs(s.cation.z) * abs(s.anion.z) * √I) / (√I_one + √I) +
    (abs(s.cation.z) * abs(s.anion.z)) / (abs(s.cation.z) + abs(s.anion.z)) * (
        F₁(s, active_salts, sss, I, W) / abs(s.cation.z) +
        F₂(s, active_salts, sss, I, W) / abs(s.anion.z)
    )
end

# Equation 7
function F₁(s::Salt, active_salts, sss, I, W)
    sum([Y(ss, I, W) * logγ⁰₁₂(ss, I) * √I_one +
         Aᵧ * √I / (√I_one + √I) * abs(ss.cation.z) * abs(ss.anion.z) * Y(ss, I, W)
         for
         ss in same_cation(s, active_salts, sss)])
end

# Equation 8
function F₂(s::Salt, active_salts, sss, I, W)
    sum([X(ss, I, W) * logγ⁰₁₂(ss, I) * √I_one +
         Aᵧ * √I / (√I_one + √I) * abs(ss.cation.z) * abs(ss.anion.z) * X(ss, I, W)
         for
         ss in same_anion(s, active_salts, sss)])
end

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

@constants T₀₂=273.15 [unit = u"K", description = "Standard temperature 2"]
@constants c_1=0.005 [
    unit = u"K^-1",
    description = "Constant for Fountoukis and Nenes (2007) Eq. 14"
]

# Equation 14
A(I) = -((0.41√I / (√I_one + √I)) + 0.039(I / I_one)^0.92)
function logγ₁₂(s::Salt, active_salts, sss, I, W)
    (1.125 - c_1 * (T - T₀₂)) * logγ₁₂T⁰(s, active_salts, sss, I, W) / √I_one -
    (0.125 - c_1 * (T - T₀₂)) * A(I)
end

"""
Calculate the activity coefficient of a salt based on Section 2.2
in Fountoukis and Nenes (2007).
"""
function activity(s::Salt, active_salts, sss, I, W)
    (s.cation.m / W)^s.ν_cation *
    (s.anion.m / W)^s.ν_anion *
    exp(logγ₁₂(s, active_salts, sss, I, W))^(s.ν_cation + s.ν_anion)
end

### Special cases

abstract type SpecialSalt <: SaltLike end

"""
The activity of a SpecialSalt is the same as for a salt except that it has
a special activity coefficient as defined in the footnotes to Table 4.
"""
function activity(s::SpecialSalt, active_salts, sss, I, W)
    (s.cation.m / W)^s.ν_cation *
    (s.anion.m / W)^s.ν_anion *
    γ₁₂(s, active_salts, sss, I, W)^(s.ν_cation + s.ν_anion)
end

logγ⁰₁₂(s::SpecialSalt, I) = abs(s.cation.z) * abs(s.anion.z) * log(Γ⁰(0, I) / I_one) # Assuming q = 0 for SpecialSalts

# Generate special salt structs for multiple dispatch.
specialsaltnames = [:CaSO4, :KHSO4, :NH4HSO4, :NaHSO4, :NH43HSO42]
for (i, name) in enumerate(specialsaltnames)
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
            drh::Any
            "Enthalpy term (-18/1000R L_s m_s)"
            l_term::Any
            "Kusik-Meissner Binary activity parameter"
            q::Number
        end
    end)
end

@mtkmodel CaSO4aqs begin
    @description "CaSO4 salt. From Table 4 footnote a, CaSO4 has an activity coefficient of zero."
    @extend Salt(; cation, ν_c, ν_a, drh, l_t, q)
end

"""
"""
γ₁₂(s::CaSO4aqs, active_salts, sss, I, W) = 1.0e-20

@mtkmodel KHSO4aqs begin
    @description "KHSO4 Salt. From Table 4 footnote b, KHSO4 has a unique activity coefficient"
    @extend Salt(; cation, ν_c, ν_a, drh, l_t, q)
end

"""
From Table 4 footnote b, KHSO4 has a unique activity coefficient
"""
function γ₁₂(s::KHSO4aqs, active_salts, sss, I, W)
    (
        exp(logγ₁₂(sss[:HHSO4], active_salts, sss, I, W)) *
        exp(logγ₁₂(sss[:KCl], active_salts, sss, I, W)) /
        exp(logγ₁₂(sss[:HCl], active_salts, sss, I, W))
    )^(1 / 2)
end

@mtkmodel NH4HSO4aqs begin
    @description "NH4HSO4 Salt. From Table 4 footnote c, NH4HSO4 has a unique activity coefficient"
    @extend Salt(; cation, ν_c, ν_a, drh, l_t, q)
end

"""
From Table 4 footnote c, NH4HSO4 has a unique activity coefficient
"""
function γ₁₂(s::NH4HSO4aqs, active_salts, sss, I, W)
    (
        exp(logγ₁₂(sss[:HHSO4], active_salts, sss, I, W)) *
        exp(logγ₁₂(sss[:NH4Cl], active_salts, sss, I, W)) /
        exp(logγ₁₂(sss[:HCl], active_salts, sss, I, W))
    )^(1 / 2)
end

@mtkmodel NaHSO4aqs begin
    @description "NaHSO4 Salt. From Table 4 footnote d, NaHSO4 has a unique activity coefficient"
    @extend Salt(; cation, ν_c, ν_a, drh, l_t, q)
end

"""
From Table 4 footnote d, NaHSO4 has a unique activity coefficient
"""
function γ₁₂(s::NaHSO4aqs, active_salts, sss, I, W)
    (
        exp(logγ₁₂(sss[:HHSO4], active_salts, sss, I, W)) *
        exp(logγ₁₂(sss[:NaCl], active_salts, sss, I, W)) /
        exp(logγ₁₂(sss[:HCl], active_salts, sss, I, W))
    )^(1 / 2)
end

@mtkmodel NH43HSO42aqs begin
    @description "NH43HSO42 Salt. From Table 4 footnote e, NH43HSO42 has a unique activity coefficient"
    @extend Salt(; cation, ν_c, ν_a, drh, l_t, q)
end

"""
From Table 4 footnote e, NH43HSO42 has a unique activity coefficient
"""
function γ₁₂(s::NH43HSO42aqs, active_salts, sss, I, W)
    (
        exp(logγ₁₂(sss[:NH42SO4], active_salts, sss, I, W))^3 *
        γ₁₂(sss[:NH4HSO4], active_salts, sss, I, W)
    )^(1 / 5)
end

"""
Create a system of equations for the ionic strength of the multicomponent solution, as specified by
Fountoukis and Nenes (2007), between equations 8 and 9: ``I = \\frac{1}{2} \\sum_i m_i z_i^2``
"""
function IonicStrength(t, all_Ions, W)
    @variables I(t)=1.0e-4 [unit = u"mol/kg_water", description = "Ionic strength"]
    I = ParentScope(I)
    @constants I_one=1 [unit = u"mol/kg_water", description = "An ionic strength of 1"]
    @named ionic_strength = ODESystem(
        [I ~ 1 / 2 * sum([ion.m / W * ion.z^2 for ion in all_Ions])], t, [I], [])
end

"""
Create an equation system for the activities of all the species listed in all_species,
where `f` is a function that takes a species as its argument and returns its activity.
"""
function Activities(t, all_species, f)
    eqs = Equation[]
    vars = []
    for s in all_species
        rhs = f(s)
        unit = ModelingToolkit.get_unit(rhs)
        x = Symbol("a_", nameof(s))
        a = only(
            @variables $x(t)=1 [unit = unit, description = "Activity of $(nameof(s))"]
        )
        push!(vars, a)
        push!(eqs, a ~ rhs)
    end
    ODESystem(eqs, t, vars, []; name = :activities)
end
