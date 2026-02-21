"An aqueous ion."
@component function Ion(; name = :Ion, logm_guess = nothing, z = 0)
    @constants begin
        _z = z, [description = "Valence (charge) of the ion (dimensionless)"]
        m_one = 1.0, [unit = u"mol/kg", description = "unit molality"]
    end
    @variables begin
        #! format: off
        z(t), [description = "Valence (charge) of the ion (dimensionless)", guess = Float64(z)]
        logm(t), [description = "Log of the molality of ion in water", guess = logm_guess]
        M(t), [description = "Molarity of ion in air", unit = u"mol/m^3", guess = 1.0e-10]
        W(t), [description = "Aerosol water content in air", unit = u"kg/m^3", guess = 1.0e-7]
        #! format: on
    end
    eqs = [
        z ~ _z,
        M ~ exp(logm) * m_one * W,
    ]
    return System(eqs, t; name)
end

"""
An aqueous salt comprised of a cation, an anion, and an activity parameter (q).
q values are given in Table 4 of Fountoukis and Nenes (2007).
"""
@component function Salt(; name = :Salt,
    is_KHSO4 = false,
    is_NH4HSO4 = false,
    is_NaHSO4 = false,
    is_NH43HSO42 = false,
    salt1 = nothing,
    salt2 = nothing,
    salt3 = nothing,
    can_precipitate = false,
    logm_guess = nothing,
    M_salt_guess = nothing,
    drh = 0.0,
    l_t = NaN,
    q = 0.0,
    ν_cation = 1,
    ν_anion = 1,
    z_cation = 1,
    z_anion = 1,
)
    # Precompute derived constants from keyword arguments
    _zz = z_cation * z_anion
    _B = 0.75 - 0.065 * q

    @constants begin
        # NOTE: The paper (between equations 6 and 7) says that the units of Aᵧ are kg^0.5 mol^−0.5, but the equation
        # doesn't work unless those units are inverted, and Bromley (1973) agrees with that.
        Aᵧ = 0.511,
            [unit = u"mol^0.5/kg^0.5", description = "Debye-Hückel constant at 298.15 K"]
        ln10 = 2.302585093,
            [description = "Natural log of 10, for converting log₁₀ constants to ln (dimensionless)"]
        T₀₂ = 273.15, [unit = u"K", description = "Standard temperature 2"]
        c_1 = 0.005,
            [unit = u"K^-1", description = "Constant for Fountoukis and Nenes (2007) Eq. 14"]

        # Salt properties
        drh = drh, [description = "Deliquescence relative humidity at 298.15K (dimensionless)"]
        l_t = l_t, [description = "Enthalpy term (-18/1000R L_s m_s) (dimensionless, or K when applicable)"]
        q = q, [description = "Kusik-Meissner Binary activity parameter (dimensionless)"]
        ν_cation = ν_cation, [description = "Number of moles of cation per mole of salt (dimensionless)"]
        ν_anion = ν_anion, [description = "Number of moles of anion per mole of salt (dimensionless)"]
        z_cation = z_cation, [description = "Absolute value of the cation charge (dimensionless)"]
        z_anion = z_anion, [description = "Absolute value of the anion charge (dimensionless)"]

        # Derived constants (precomputed from keyword arguments)
        B = _B, [description = "Parameter for Kusik-Meissner binary activity (dimensionless)"]
        zz = _zz, [description = "Product of the absolute values of the cation and anion charges (dimensionless)"]

        # Unit conversions
        I_one = 1, [unit = u"mol/kg", description = "An ionic strength of 1"]
        m_one = 1, [unit = u"mol/kg", description = "A molality of 1"]
        M_zero = 0, [unit = u"mol/m^3", description = "A molarity of 0"]
    end
    @variables begin
        #! format: off
        T(t), [description = "Temperature", unit = u"K", guess = 293.15]
        RH(t), [description = "Relative humidity (0-1)", guess = 0.3]
        logm(t), [description = "log of the molality of the salt in water", guess = logm_guess]
        M(t), [description = "Aqueous molarity of the salt in air", unit = u"mol/m^3", guess = M_salt_guess === nothing ? 1.0e-10 : M_salt_guess]
        M_precip(t), [description = "Precipitated solid", unit = u"mol/m^3", guess = 0.0]
        deliquesced(t), [description = "Whether the salt is deliquesced (1) or not (0) (dimensionless)", guess = 0.5]
        W(t), [description = "Aerosol water content in air", unit = u"kg/m^3", guess = 1.0e-7]
        X(t), [description = "Parameter for activity coefficient calculation (dimensionless)", guess = 0.5]
        Y(t), [description = "Parameter for activity coefficient calculation (dimensionless)", guess = 0.5]
        I(t), [description = "Ionic strength", unit = u"mol/kg", guess = 1]
        logγ⁰(t), [description = "Log of the mean ionic activity coefficient for the single solute solution (dimensionless)", guess = 0.0]
        logγₜ₀(t), [description = "Log of the multi-component activity coefficient at 298.15K (dimensionless)", guess = 0]
        logγ(t), [description = "Log of the activity coefficient (dimensionless)", guess = 0]
        logΓ⁰(t), [description = "Parameter for Kusik-Meissner activity coefficient (dimensionless)", guess = 0.0]
        logΓ⁺(t), [description = "Parameter for Pitzer activity coefficient (dimensionless)", guess = 0.0]
        C(t), [description = "Parameter for Pitzer activity coefficient (dimensionless)", guess = 1.0]
        F_cat(t), [description = "Activity contribution from the cation (dimensionless)", guess = 0.0]
        F_an(t), [description = "Activity contribution from the anion (dimensionless)", guess = 0.0]
        A(t), [description = "Parameter used in Equation 14 (dimensionless)", guess = 0.0]
        loga(t), [description = "Log of the activity of the salt (dimensionless)", guess = 0.0]
        #! format: on
    end

    if can_precipitate
        @variables begin
            M₀(t), [description = "Candidate molarity (internal variable)", unit = u"mol/m^3", guess = 1.0e-10]
        end
    end

    eqs = Equation[
        # Equation 6 (constants multiplied by ln10 to convert from log₁₀ to natural log)
        logγₜ₀ ~
            -ln10 * Aᵧ * zz * √I / (√I_one + √I) / √I_one + # NOTE: The last √I_one here is not in the paper but is needed to make the units balance.
            (zz / (z_cation + z_anion)) * (F_cat / z_cation + F_an / z_anion),

        # Supplemental equations after equations 7 and 8
        Y ~ ((z_cation + z_anion) / 2)^2 * ν_anion * exp(logm) * m_one / I,
        X ~ ((z_cation + z_anion) / 2)^2 * ν_cation * exp(logm) * m_one / I,
        # Equation 9
        logγ⁰ ~ zz * logΓ⁰,
        # Equation 10
        logΓ⁰ ~ log(1 + B * (1 + 0.1I / I_one)^q - B) + logΓ⁺,
        # Equation 12 (constant multiplied by ln10 to convert from log₁₀ to natural log)
        logΓ⁺ ~ -0.5107 * ln10 * √I / (√I_one + C * √I),
        # Equation 13
        C ~ 1 + 0.055q * exp(-0.023I^3 / I_one^3),

        # Equation in text below Equation 14 (constants multiplied by ln10 to convert from log₁₀ to natural log)
        A ~ -ln10 * ((0.41√I / (√I_one + √I)) + 0.039(I / I_one)^0.92),

        # Activity (Section 2.2)
        # ln(a) = ln(ν_c^ν_c * ν_a^ν_a) + ν*(ln(m) + ln(γ±))
        logm ~
            (loga - log(ν_cation^ν_cation * ν_anion^ν_anion) - (ν_cation + ν_anion) * logγ) /
            (ν_cation + ν_anion),

        M_precip ~ 0,
        deliquesced ~ (tanh((RH - drh) * 30) + 1) / 2,
    ]

    # Equation 14 (conditional)
    if is_KHSO4 || is_NH4HSO4 || is_NaHSO4  # From Table 4 footnote b, c & d
        push!(eqs,
            logγ ~
                0.5 *
                (ParentScope(salt1.logγ) + ParentScope(salt2.logγ) - ParentScope(salt3.logγ)))
    elseif is_NH43HSO42
        push!(eqs,
            logγ ~ (ParentScope(salt1.logγ) * 3 + ParentScope(salt2.logγ)) * (1 / 5)) # From Table 4 footnote e
    else
        push!(eqs,
            logγ ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀ - (0.125 - c_1 * (T - T₀₂)) * A)
    end

    # M equation (conditional on can_precipitate)
    if can_precipitate
        push!(eqs,
            M ~ (1 - deliquesced) * min(
                M₀, # Maximum dissolved concentration
                exp(logm) * m_one * W, # Equilibrium with solid.
            ) + deliquesced * M₀) # Fully dissolved above the DRH
    else
        push!(eqs, M ~ exp(logm) * m_one * W)
    end

    return System(eqs, t; name)
end

"""Water content for a binary aqueous aerosol solution using data from
Fountoukis and Nenes (2007) Table 7, for use in Equation 16."""
@component function BinaryMolality(; name = :BinaryMolality, k_0, k_1, k_2, k_3, k_4, k_5, k_6)
    @constants begin
        k_0 = k_0, [description = "polynomial fit coefficient (dimensionless)"]
        k_1 = k_1, [description = "polynomial fit coefficient (dimensionless)"]
        k_2 = k_2, [description = "polynomial fit coefficient (dimensionless)"]
        k_3 = k_3, [description = "polynomial fit coefficient (dimensionless)"]
        k_4 = k_4, [description = "polynomial fit coefficient (dimensionless)"]
        k_5 = k_5, [description = "polynomial fit coefficient (dimensionless)"]
        k_6 = k_6, [description = "polynomial fit coefficient (dimensionless)"]
        m_one = 1.0, [unit = u"mol/kg", description = "unit molality"]
    end
    @variables begin
        RH(t), [description = "Relative humidity (dimensionless)", guess = 0.3]
        m_aw(t), [description = "molality of the binary solution", unit = u"mol/kg", guess = 1.0]
    end
    eqs = [
        m_aw ~ (k_0 + k_1 * RH + k_2 * RH^2 + k_3 * RH^3 + k_4 * RH^4 + k_5 * RH^5) * m_one,
    ]
    return System(eqs, t; name)
end

"Aqueous behavior"
@component function Aqueous(; name = :Aqueous, q0 = 0.0)
    @constants begin
        # NOTE: The paper (between equations 6 and 7) says that the units of Aᵧ are kg^0.5 mol^−0.5, but the equation
        # doesn't work unless those units are inverted, and Bromley (1973) agrees with that.
        Aᵧ = 0.511,
            [unit = u"mol^0.5/kg^0.5", description = "Debye-Hückel constant at 298.15 K"]
        ln10 = 2.302585093,
            [description = "Natural log of 10, for converting log₁₀ constants to ln (dimensionless)"]
        I_one = 1, [unit = u"mol/kg", description = "An ionic strength of 1"]

        m_one = 1.0, [unit = u"mol/kg", description = "unit molality"]
        w_c = 1.0e-6, [unit = u"kg/m^3", description = "Water content in air"]
        W_min = 1.0e-7, [unit = u"kg/m^3", description = "Minimum water content in air"]
    end
    @variables begin
        T(t), [description = "Temperature", unit = u"K", guess = 293.15]
        RH(t), [description = "Relative humidity (0-1)", guess = 0.3]
        I(t), [description = "Ionic strength", unit = u"mol/kg", guess = 1]
        W(t),
            [description = "Aerosol water content in air", unit = u"kg/m^3", guess = 1.0e-7]
        W_x(t),
            [description = "Aerosol water content in air", unit = u"kg/m^3", guess = 1.0e-7]

        F_Ca(t), [description = "Activity contribution from Ca-containing salts", guess = 0.0]
        F_K(t), [description = "Activity contribution from K-containing salts", guess = 0.0]
        F_Mg(t), [description = "Activity contribution from Mg-containing salts", guess = 0.0]
        F_Na(t), [description = "Activity contribution from Na-containing salts", guess = 0.0]
        F_NH4(t), [description = "Activity contribution from NH4-containing salts", guess = 0.0]
        F_Cl(t), [description = "Activity contribution from Cl-containing salts", guess = 0.0]
        F_H(t), [description = "Activity contribution from H-containing salts", guess = 0.0]
        F_NO3(t), [description = "Activity contribution from NO3-containing salts", guess = 0.0]
        F_SO4(t), [description = "Activity contribution from SO4-containing salts", guess = 0.0]
        F_HSO4(t), [description = "Activity contribution from HSO4-containing salts", guess = 0.0]
        F_OH(t), [description = "Activity contribution from OH-containing salts", guess = 0.0]

        Aᵧ_term(t), [description = "Debye-Hückel term used in Equation 7 and 8 (dimensionless)", guess = 0.0]
    end

    # Cations
    H = Ion(; name = :H, z = 1, logm_guess = 0)
    # Anions
    OH = Ion(; name = :OH, z = abs(-1), logm_guess = 0)
    # Neutral species
    NH3 = Ion(; name = :NH3, z = 0, logm_guess = 0)
    HNO3_aq = Ion(; name = :HNO3_aq, z = 0, logm_guess = 0)
    HCl_aq = Ion(; name = :HCl_aq, z = 0, logm_guess = 0)

    # Salts
    CaNO32 = Salt(; name = :CaNO32,
        z_cation = 2, ν_anion = 2, drh = 0.4906, l_t = 509.4, q = 0.93,
        logm_guess = 0, M_salt_guess = 1.9e-8, can_precipitate = false)
    CaCl2 = Salt(; name = :CaCl2,
        z_cation = 2, ν_anion = 2, drh = 0.283, l_t = 551.1, q = 2.4,
        logm_guess = 0, can_precipitate = false)
    CaSO4 = Salt(; name = :CaSO4,
        z_cation = 2, z_anion = 2, drh = 0.97, l_t = NaN, q = q0,
        logm_guess = 0, can_precipitate = false)
    K2SO4 = Salt(; name = :K2SO4,
        ν_cation = 2, z_anion = 2, drh = 0.9751, l_t = 35.6, q = -0.25,
        logm_guess = 0, can_precipitate = true)
    KNO3 = Salt(; name = :KNO3,
        z_cation = 1, z_anion = 1, drh = 0.9248, l_t = NaN, q = -2.33,
        logm_guess = 0)
    KCl = Salt(; name = :KCl,
        z_cation = 1, z_anion = 1, drh = 0.8426, l_t = 158.9, q = 0.92,
        logm_guess = 0)
    MgSO4 = Salt(; name = :MgSO4,
        z_cation = 2, z_anion = 2, drh = 0.8613, l_t = -714.5, q = 0.15,
        logm_guess = 0)
    MgNO32 = Salt(; name = :MgNO32,
        z_cation = 2, ν_anion = 2, drh = 0.54, l_t = 230.2, q = 2.32,
        logm_guess = 0, can_precipitate = false)
    MgCl2 = Salt(; name = :MgCl2,
        z_cation = 2, ν_anion = 2, drh = 0.3284, l_t = 42.23, q = 2.9,
        logm_guess = 0, can_precipitate = false)
    NaCl = Salt(; name = :NaCl,
        z_cation = 1, z_anion = 1, drh = 0.7528, l_t = 25.0, q = 2.23,
        logm_guess = 0)
    Na2SO4 = Salt(; name = :Na2SO4,
        ν_cation = 2, z_anion = 2, drh = 0.93, l_t = 80.0, q = -0.19,
        logm_guess = 0, can_precipitate = true)
    NaNO3 = Salt(; name = :NaNO3,
        z_cation = 1, z_anion = 1, drh = 0.7379, l_t = 304.0, q = -0.39,
        logm_guess = 0)
    NH42SO4 = Salt(; name = :NH42SO4,
        ν_cation = 2, z_anion = 2, drh = 0.7997, l_t = 80.0, q = -0.25,
        logm_guess = 0)
    NH4NO3 = Salt(; name = :NH4NO3,
        z_cation = 1, z_anion = 1, drh = 0.6183, l_t = 852.0, q = -1.15,
        logm_guess = 0)
    NH4Cl = Salt(; name = :NH4Cl,
        z_cation = 1, z_anion = 1, drh = 0.771, l_t = 239.0, q = 0.82,
        logm_guess = 0)
    HHSO4 = Salt(; name = :HHSO4,
        z_cation = 1, z_anion = 1, drh = 0.0, l_t = NaN, q = 8.0)
    H2SO4 = Salt(; name = :H2SO4,
        z_cation = 1, z_anion = 1, drh = 0.0, l_t = NaN, q = -0.1,
        logm_guess = 0)
    HNO3 = Salt(; name = :HNO3,
        z_cation = 1, z_anion = 1, drh = 0, l_t = NaN, q = 2.6,
        logm_guess = 0)
    HCl = Salt(; name = :HCl,
        z_cation = 1, z_anion = 1, drh = 0, l_t = NaN, q = 6.0,
        logm_guess = 0)
    KHSO4 = Salt(; name = :KHSO4,
        z_cation = 1, z_anion = 1, drh = 0.86, l_t = NaN, q = q0,
        is_KHSO4 = true, salt1 = HHSO4, salt2 = KCl, salt3 = HCl,
        logm_guess = 0)
    NH4HSO4 = Salt(; name = :NH4HSO4,
        z_cation = 1, z_anion = 1, drh = 0.4, l_t = 384.0, q = q0,
        is_NH4HSO4 = true, salt1 = HHSO4, salt2 = NH4Cl, salt3 = HCl,
        logm_guess = 0)
    NaHSO4 = Salt(; name = :NaHSO4,
        z_cation = 1, z_anion = 1, drh = 0.52, l_t = -45.0, q = q0,
        is_NaHSO4 = true, salt1 = HHSO4, salt2 = NaCl, salt3 = HCl,
        logm_guess = 0)
    NH43HSO42 = Salt(; name = :NH43HSO42,
        ν_cation = 3, ν_anion = 2, drh = 0.69, l_t = 186.0, q = q0,
        is_NH43HSO42 = true, salt1 = NH42SO4, salt2 = NH4HSO4,
        logm_guess = 0)

    # Species that are not in the paper. Assume q=0
    HSO4_dissociated = Salt(; name = :HSO4_dissociated,
        z_anion = abs(-2), drh = 0, l_t = NaN, q = 0.0, logm_guess = 1) # H + SO4
    NH3_dissociated = Salt(; name = :NH3_dissociated,
        drh = 0, l_t = NaN, q = 0.0, logm_guess = 1) # NH4 + OH
    H2O_dissociated = Salt(; name = :H2O_dissociated,
        drh = 0, l_t = NaN, q = 0.0, logm_guess = 1) # H + OH

    # Water content
    #! format: off
    maw_CaNO32 = BinaryMolality(; name = :maw_CaNO32,
        k_0 = 36.356, k_1 = -165.66, k_2 = 447.46, k_3 = -673.55, k_4 = 510.91, k_5 = -155.56, k_6 = 0)
    maw_CaCl2 = BinaryMolality(; name = :maw_CaCl2,
        k_0 = 20.847, k_1 = -97.599, k_2 = 273.220, k_3 = -422.120, k_4 = 331.160, k_5 = -105.450, k_6 = 0)
    maw_KHSO4 = BinaryMolality(; name = :maw_KHSO4,
        k_0 = 1.061, k_1 = -0.101, k_2 = 1.579e-2, k_3 = -1.950e-3, k_4 = 9.515e-5, k_5 = -1.547e-6, k_6 = 0)
    maw_K2SO4 = BinaryMolality(; name = :maw_K2SO4,
        k_0 = 1061.51, k_1 = -4748.97, k_2 = 8096.16, k_3 = -6166.16, k_4 = 1757.47, k_5 = 0, k_6 = 0)
    maw_KNO3 = BinaryMolality(; name = :maw_KNO3,
        k_0 = 1.2141e4, k_1 = -5.1173e4, k_2 = 8.1252e4, k_3 = -5.7527e4, k_4 = 1.5305e4, k_5 = 0, k_6 = 0)
    maw_KCl = BinaryMolality(; name = :maw_KCl,
        k_0 = 179.721, k_1 = -721.266, k_2 = 1161.03, k_3 = -841.479, k_4 = 221.943, k_5 = 0, k_6 = 0)
    maw_MgSO4 = BinaryMolality(; name = :maw_MgSO4,
        k_0 = -0.778, k_1 = 177.74, k_2 = -719.79, k_3 = 1174.6, k_4 = -863.44, k_5 = 232.31, k_6 = 0)
    maw_MgNO32 = BinaryMolality(; name = :maw_MgNO32,
        k_0 = 12.166, k_1 = -16.154, k_2 = 0, k_3 = 10.886, k_4 = 0, k_5 = -6.815, k_6 = 0)
    maw_MgCl2 = BinaryMolality(; name = :maw_MgCl2,
        k_0 = 11.505, k_1 = -26.518, k_2 = 34.937, k_3 = -19.829, k_4 = 0, k_5 = 0, k_6 = 0)
    maw_NaNO3 = BinaryMolality(; name = :maw_NaNO3,
        k_0 = 0.9988, k_1 = -2.6947e-2, k_2 = 1.9610e-4, k_3 = 2.8154e-5, k_4 = 6.1359e-7, k_5 = 0, k_6 = 0)
    maw_NaHSO4 = BinaryMolality(; name = :maw_NaHSO4,
        k_0 = 1.0614, k_1 = -0.1014, k_2 = 1.5796e-2, k_3 = -1.9501e-3, k_4 = 9.5147e-5, k_5 = -1.5473e-6, k_6 = 0)
    maw_NaCl = BinaryMolality(; name = :maw_NaCl,
        k_0 = 1.0084, k_1 = -4.9390e-2, k_2 = 8.888e-3, k_3 = -2.1570e-3, k_4 = 1.6170e-4, k_5 = 1.99e-6, k_6 = -1.142e-7)
    maw_Na2SO4 = BinaryMolality(; name = :maw_Na2SO4,
        k_0 = 1.0052, k_1 = -6.4840e-2, k_2 = 3.519e-2, k_3 = -1.3190e-2, k_4 = 1.9250e-3, k_5 = -1.224e-4, k_6 = 2.87e-6)
    maw_NH42SO4 = BinaryMolality(; name = :maw_NH42SO4,
        k_0 = 0.9968, k_1 = -2.9690e-2, k_2 = 1.735e-5, k_3 = -3.2530e-4, k_4 = 3.5710e-5, k_5 = -9.7870e-7, k_6 = 0)
    maw_NH4Cl = BinaryMolality(; name = :maw_NH4Cl,
        k_0 = 0.9968, k_1 = -2.6110e-2, k_2 = -1.5990e-3, k_3 = 1.3550e-4, k_4 = -2.317e-6, k_5 = -1.113e-8, k_6 = 0)
    maw_NH4NO3 = BinaryMolality(; name = :maw_NH4NO3,
        k_0 = 1.0053, k_1 = -2.4991e-2, k_2 = 4.4688e-4, k_3 = 1.6453e-5, k_4 = -3.8940e-7, k_5 = -4.7668e-8, k_6 = 1.3753e-9)
    maw_NH4HSO4 = BinaryMolality(; name = :maw_NH4HSO4,
        k_0 = 1.0261, k_1 = -4.9766e-2, k_2 = 3.2757e-3, k_3 = -2.4477e-4, k_4 = 1.0766e-5, k_5 = -1.8329e-7, k_6 = 0)
    maw_NH43HSO42 = BinaryMolality(; name = :maw_NH43HSO42,
        k_0 = 1.0088, k_1 = -5.3730e-2, k_2 = 1.4201e-3, k_3 = -9.2484e-4, k_4 = 2.2796e-4, k_5 = -1.5445e-5, k_6 = 0)
    #! format: on

    # Collect all subsystems
    all_salts = [CaNO32, CaCl2, CaSO4, K2SO4, KNO3, KCl, MgSO4, MgNO32, MgCl2,
        NaCl, Na2SO4, NaNO3, NH42SO4, NH4NO3, NH4Cl, HHSO4, H2SO4, HNO3, HCl,
        KHSO4, NH4HSO4, NaHSO4, NH43HSO42, HSO4_dissociated, NH3_dissociated, H2O_dissociated]
    all_ions = [H, OH, NH3, HNO3_aq, HCl_aq]
    all_maws = [maw_CaNO32, maw_CaCl2, maw_KHSO4, maw_K2SO4, maw_KNO3, maw_KCl,
        maw_MgSO4, maw_MgNO32, maw_MgCl2, maw_NaNO3, maw_NaHSO4, maw_NaCl,
        maw_Na2SO4, maw_NH42SO4, maw_NH4Cl, maw_NH4NO3, maw_NH4HSO4, maw_NH43HSO42]

    eqs = [
        Aᵧ_term ~ ln10 * Aᵧ * √I / (√I_one + √I) / √I_one, # NOTE: The last √I_one here is not in the paper but is needed to make the units balance. Constants multiplied by ln10 to convert from log₁₀ to natural log.
        # Equation 7
        F_NH4 ~
            sum([s.Y * s.logγ⁰ for s in [NH4NO3, NH4Cl, NH4HSO4, NH42SO4, NH43HSO42]]) + Aᵧ_term * sum([s.zz * s.Y for s in [NH4NO3, NH4Cl, NH4HSO4, NH42SO4, NH43HSO42]]),
        F_Na ~
            sum([s.Y * s.logγ⁰ for s in [NaCl, Na2SO4, NaNO3, NaHSO4]]) + Aᵧ_term * sum([s.zz * s.Y for s in [NaCl, Na2SO4, NaNO3, NaHSO4]]),
        F_H ~
            sum([s.Y * s.logγ⁰ for s in [HCl, HNO3]]) + Aᵧ_term * sum([s.zz * s.Y for s in [HCl, HNO3]]),
        F_Ca ~
            sum([s.Y * s.logγ⁰ for s in [CaNO32, CaCl2, CaSO4]]) + Aᵧ_term * sum([s.zz * s.Y for s in [CaNO32, CaCl2, CaSO4]]),
        F_K ~
            sum([s.Y * s.logγ⁰ for s in [KHSO4, K2SO4, KNO3, KCl]]) + Aᵧ_term * sum([s.zz * s.Y for s in [KHSO4, K2SO4, KNO3, KCl]]),
        F_Mg ~
            sum([s.Y * s.logγ⁰ for s in [MgSO4, MgNO32, MgCl2]]) + Aᵧ_term * sum([s.zz * s.Y for s in [MgSO4, MgNO32, MgCl2]]),

        # Equation 8
        F_Cl ~
            sum([s.X * s.logγ⁰ for s in [NaCl, KCl, MgCl2, CaCl2, NH4Cl, HCl]]) + Aᵧ_term * sum([s.zz * s.X for s in [NaCl, KCl, MgCl2, CaCl2, NH4Cl, HCl]]),
        F_NO3 ~
            sum([s.X * s.logγ⁰ for s in [NaNO3, KNO3, MgNO32, CaNO32, NH4NO3, HNO3]]) + Aᵧ_term * sum([s.zz * s.X for s in [NaNO3, KNO3, MgNO32, CaNO32, NH4NO3, HNO3]]),
        F_SO4 ~
            sum([s.X * s.logγ⁰ for s in [Na2SO4, K2SO4, MgSO4, CaSO4, NH42SO4]]) + Aᵧ_term * sum([s.zz * s.X for s in [Na2SO4, K2SO4, MgSO4, CaSO4, NH42SO4]]),
        F_HSO4 ~
            sum([s.X * s.logγ⁰ for s in [KHSO4, NaHSO4, NH4HSO4, NH43HSO42]]) + Aᵧ_term * sum([s.zz * s.X for s in [KHSO4, NaHSO4, NH4HSO4, NH43HSO42]]),
        F_OH ~
            sum([s.X * s.logγ⁰ for s in [H2O_dissociated, NH3_dissociated]]) + Aᵧ_term * sum([s.zz * s.X for s in [H2O_dissociated, NH3_dissociated]]),

        CaNO32.F_cat ~ F_Ca,
        CaNO32.F_an ~ F_NO3,
        CaCl2.F_cat ~ F_Ca,
        CaCl2.F_an ~ F_Cl,
        CaSO4.F_cat ~ F_Ca,
        CaSO4.F_an ~ F_SO4,
        KHSO4.F_cat ~ F_K,
        KHSO4.F_an ~ F_HSO4,
        K2SO4.F_cat ~ F_K,
        K2SO4.F_an ~ F_SO4,
        KNO3.F_cat ~ F_K,
        KNO3.F_an ~ F_NO3,
        KCl.F_cat ~ F_K,
        KCl.F_an ~ F_Cl,
        MgSO4.F_cat ~ F_Mg,
        MgSO4.F_an ~ F_SO4,
        MgNO32.F_cat ~ F_Mg,
        MgNO32.F_an ~ F_NO3,
        MgCl2.F_cat ~ F_Mg,
        MgCl2.F_an ~ F_Cl,
        NaCl.F_cat ~ F_Na,
        NaCl.F_an ~ F_Cl,
        Na2SO4.F_cat ~ F_Na,
        Na2SO4.F_an ~ F_SO4,
        NaNO3.F_cat ~ F_Na,
        NaNO3.F_an ~ F_NO3,
        NH42SO4.F_cat ~ F_NH4,
        NH42SO4.F_an ~ F_SO4,
        NH4NO3.F_cat ~ F_NH4,
        NH4NO3.F_an ~ F_NO3,
        NH4Cl.F_cat ~ F_NH4,
        NH4Cl.F_an ~ F_Cl,
        NH4HSO4.F_cat ~ F_NH4,
        NH4HSO4.F_an ~ F_HSO4,
        NaHSO4.F_cat ~ F_Na,
        NaHSO4.F_an ~ F_HSO4,
        NH43HSO42.F_cat ~ F_NH4,
        NH43HSO42.F_an ~ F_HSO4,
        HHSO4.F_cat ~ F_H,
        HHSO4.F_an ~ F_HSO4,
        H2SO4.F_cat ~ F_H,
        H2SO4.F_an ~ F_HSO4,
        HNO3.F_cat ~ F_H,
        HNO3.F_an ~ F_NO3,
        HCl.F_cat ~ F_H,
        HCl.F_an ~ F_Cl,
        NH3_dissociated.F_cat ~ F_NH4,
        NH3_dissociated.F_an ~ F_OH,
        H2O_dissociated.F_cat ~ F_H,
        H2O_dissociated.F_an ~ F_OH,
        HSO4_dissociated.F_cat ~ F_H,
        HSO4_dissociated.F_an ~ F_SO4,

        # Ionic strength (equation in text below Equation 8)
        I ~
            0.5 *
            (
            sum(
                [
                    exp(salt.logm) * salt.ν_cation * salt.z_cation^2 +
                        exp(salt.logm) * salt.ν_anion * salt.z_anion^2 for salt in [
                            CaNO32,
                            CaCl2,
                            CaSO4,
                            KHSO4,
                            K2SO4,
                            KNO3,
                            KCl,
                            MgSO4,
                            MgNO32,
                            MgCl2,
                            NaCl,
                            Na2SO4,
                            NaNO3,
                            NH42SO4,
                            NH4NO3,
                            NH4Cl,
                            NH4HSO4,
                            NaHSO4,
                            NH43HSO42,
                        ]
                ]
            ) +
                exp(H.logm) * H.z^2 +
                exp(OH.logm) * OH.z^2
        ) *
            m_one,

        CaNO32.I ~ I,
        CaCl2.I ~ I,
        CaSO4.I ~ I,
        KHSO4.I ~ I,
        K2SO4.I ~ I,
        KNO3.I ~ I,
        KCl.I ~ I,
        MgSO4.I ~ I,
        MgNO32.I ~ I,
        MgCl2.I ~ I,
        NaCl.I ~ I,
        Na2SO4.I ~ I,
        NaNO3.I ~ I,
        NH42SO4.I ~ I,
        NH4NO3.I ~ I,
        NH4Cl.I ~ I,
        NH4HSO4.I ~ I,
        NaHSO4.I ~ I,
        NH43HSO42.I ~ I,
        HHSO4.I ~ I,
        H2SO4.I ~ I,
        HNO3.I ~ I,
        HCl.I ~ I,
        NH3_dissociated.I ~ I,
        H2O_dissociated.I ~ I,
        HSO4_dissociated.I ~ I,

        # Water content (Section 2.3)
        W ~ max(
            CaNO32.M / maw_CaNO32.m_aw +
                CaCl2.M / maw_CaCl2.m_aw +
                KHSO4.M / maw_KHSO4.m_aw +
                K2SO4.M / maw_K2SO4.m_aw +
                KNO3.M / maw_KNO3.m_aw +
                KCl.M / maw_KCl.m_aw +
                MgSO4.M / maw_MgSO4.m_aw +
                MgNO32.M / maw_MgNO32.m_aw +
                MgCl2.M / maw_MgCl2.m_aw +
                NaCl.M / maw_NaCl.m_aw +
                Na2SO4.M / maw_Na2SO4.m_aw +
                NaNO3.M / maw_NaNO3.m_aw +
                NH42SO4.M / maw_NH42SO4.m_aw +
                NH4NO3.M / maw_NH4NO3.m_aw +
                NaHSO4.M / maw_NaHSO4.m_aw +
                NH4Cl.M / maw_NH4Cl.m_aw +
                NH4HSO4.M / maw_NH4HSO4.m_aw +
                NH43HSO42.M / maw_NH43HSO42.m_aw,
            W_min,
        ),

        maw_CaNO32.RH ~ RH,
        maw_CaCl2.RH ~ RH,
        maw_KHSO4.RH ~ RH,
        maw_K2SO4.RH ~ RH,
        maw_KNO3.RH ~ RH,
        maw_KCl.RH ~ RH,
        maw_MgSO4.RH ~ RH,
        maw_MgNO32.RH ~ RH,
        maw_MgCl2.RH ~ RH,
        maw_NaNO3.RH ~ RH,
        maw_NaHSO4.RH ~ RH,
        maw_NaCl.RH ~ RH,
        maw_Na2SO4.RH ~ RH,
        maw_NH42SO4.RH ~ RH,
        maw_NH4Cl.RH ~ RH,
        maw_NH4NO3.RH ~ RH,
        maw_NH4HSO4.RH ~ RH,
        maw_NH43HSO42.RH ~ RH,

        # Connect W to all salts and ions
        CaNO32.W ~ W,
        CaCl2.W ~ W,
        CaSO4.W ~ W,
        KHSO4.W ~ W,
        K2SO4.W ~ W,
        KNO3.W ~ W,
        KCl.W ~ W,
        MgSO4.W ~ W,
        MgNO32.W ~ W,
        MgCl2.W ~ W,
        NaCl.W ~ W,
        Na2SO4.W ~ W,
        NaNO3.W ~ W,
        NH42SO4.W ~ W,
        NH4NO3.W ~ W,
        NH4Cl.W ~ W,
        NH4HSO4.W ~ W,
        NaHSO4.W ~ W,
        NH43HSO42.W ~ W,
        HHSO4.W ~ W,
        H2SO4.W ~ W,
        HNO3.W ~ W,
        HCl.W ~ W,

        H.W ~ W,
        HSO4_dissociated.W ~ W,
        OH.W ~ W,
        NH3.W ~ W,
        NH3_dissociated.W ~ W,
        HNO3_aq.W ~ W,
        HCl_aq.W ~ W,
        H2O_dissociated.W ~ W,

        # NH4Cl and NH4NO3 precipitate directly from gas, so we assume
        # the aqueous concentration is zero.

        # HHSO4 is only used to calculate activity coefficients of other salts.
        HHSO4.logm ~ -40,

        # Mass balance
        H.M ~ sum([HCl.M, H2SO4.M, HSO4_dissociated.M, H2O_dissociated.M, HNO3.M]),
        OH.M ~ NH3_dissociated.M + H2O_dissociated.M,

        # Connect T and RH to all salts
        CaNO32.T ~ T,
        CaNO32.RH ~ RH,
        CaCl2.T ~ T,
        CaCl2.RH ~ RH,
        CaSO4.T ~ T,
        CaSO4.RH ~ RH,
        K2SO4.T ~ T,
        K2SO4.RH ~ RH,
        KNO3.T ~ T,
        KNO3.RH ~ RH,
        KCl.T ~ T,
        KCl.RH ~ RH,
        MgSO4.T ~ T,
        MgSO4.RH ~ RH,
        MgNO32.T ~ T,
        MgNO32.RH ~ RH,
        MgCl2.T ~ T,
        MgCl2.RH ~ RH,
        NaCl.T ~ T,
        NaCl.RH ~ RH,
        Na2SO4.T ~ T,
        Na2SO4.RH ~ RH,
        NaNO3.T ~ T,
        NaNO3.RH ~ RH,
        NH42SO4.T ~ T,
        NH42SO4.RH ~ RH,
        NH4NO3.T ~ T,
        NH4NO3.RH ~ RH,
        NH4Cl.T ~ T,
        NH4Cl.RH ~ RH,
        HHSO4.T ~ T,
        HHSO4.RH ~ RH,
        H2SO4.T ~ T,
        H2SO4.RH ~ RH,
        HNO3.T ~ T,
        HNO3.RH ~ RH,
        HCl.T ~ T,
        HCl.RH ~ RH,
        # Note: KHSO4, NH4HSO4, NaHSO4, NH43HSO42 use special logγ formulas
        # that don't depend on T, so T is not an unknown in those systems.
        KHSO4.RH ~ RH,
        NH4HSO4.RH ~ RH,
        NaHSO4.RH ~ RH,
        NH43HSO42.RH ~ RH,
        HSO4_dissociated.T ~ T,
        HSO4_dissociated.RH ~ RH,
        NH3_dissociated.T ~ T,
        NH3_dissociated.RH ~ RH,
        H2O_dissociated.T ~ T,
        H2O_dissociated.RH ~ RH,
    ]

    return System(eqs, t; systems = [all_salts; all_ions; all_maws], name)
end

"xxx"
