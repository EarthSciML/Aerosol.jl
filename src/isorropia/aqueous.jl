@mtkmodel Ion begin
    @description "An aqueous ion."
    @parameters begin
        z, [description = "Valence (charge) of the ion"]
    end
    @variables begin
        m(t), [description = "Molality of ion in water", unit = u"mol/kg", guess = 1]
        a(t),
        [
            description = "Activity of the ion. The activity coefficient of an ion is assumed to be one (Fountoukis and Nenes (2007), Section 3.3).",
            unit = u"mol/kg"]
    end
    @equations begin
        a ~ m
    end
end

@mtkmodel Salt begin
    @description """
An aqueous salt comprised of a cation, an anion, and an activity parameter (q).
q values are given in Table 4 of Fountoukis and Nenes (2007).
"""
    @structural_parameters begin
        cation = Ion(; name = :cation, z = 1)
        anion = Ion(; name = :anion, z = 1)
        is_CaSO4 = false # Is this salt CaSO4?
        is_KHSO4 = false # Is this salt KHSO4?
        is_NH4HSO4 = false # Is this salt NH4HSO4?
        is_NaHSO4 = false # Is this salt NaHSO4?
        is_NH43HSO42 = false # Is this salt NH43HSO42?
        salt1 = nothing
        salt2 = nothing
        salt3 = nothing
    end
    @constants begin
        # NOTE: The paper (between equations 6 and 7) says that the units of Aᵧ are kg^0.5 mol^−0.5, but the equation
        # doesn't work unless those units are inverted, and Bromley (1973) agrees with that.
        Aᵧ = 0.511,
        [
            unit = u"mol^0.5/kg^0.5",
            description = "Debye-Hückel constant at 298.15 K"
        ]
        T₀₂ = 273.15, [unit = u"K", description = "Standard temperature 2"]
        c_1 = 0.005,
        [
            unit = u"K^-1",
            description = "Constant for Fountoukis and Nenes (2007) Eq. 14"
        ]

        # Salt properties
        drh, [description = "Deliquescence relative humidity at 298.15K"]
        l_t, [description = "Enthalpy term (-18/1000R L_s m_s)"]
        q, [description = "Kusik-Meissner Binary activity parameter"]
        ν_cation = 1, [description = "Number of moles of cation per mole of salt"]
        ν_anion = 1, [description = "Number of moles of anion per mole of salt"]

        # Derived constants
        B
        zz, [description = "Product of the absolute values of the cation and anion charges"]

        # Unit conversions
        I_one = 1, [unit = u"mol/kg", description = "An ionic strength of 1"]
        m_one = 1, [unit = u"mol/kg", description = "A molality of 1"]
    end
    @parameters begin
        T = 293.15, [description = "Temperature", unit = u"K"]
    end
    @variables begin
        M(t), [description = "Molarity of the salt in air", unit = u"mol/m^3", guess = 1e-8]
        X(t)
        Y(t)
        I(t), [description = "Ionic strength", unit = u"mol/kg", guess = 1e-8]
        #! format: off
        logγ⁰(t), [description = "Log of the mean ionic activity coefficient for the single solute solution", guess=1]
        logγₜ₀(t), [description = "Log of the multi-component activity coefficient at 298.15K", guess=1]
        logγ(t), [description = "Log of the activity coefficient", guess=1]
        #! format: on
        logΓ⁰(t)
        logΓ⁺(t)
        C(t)
        F_cat(t), [description = "Activity contribution from the cation"]
        F_an(t), [description = "Activity contribution from the anion"]
        A(t), [description = "Parameter used in Equation 14"]
        loga(t), [description = "Log of the activity of the salt", guess = 0]
    end
    @equations begin
        zz ~ ParentScope(cation.z) * ParentScope(anion.z)

        # Equation 6
        logγₜ₀ ~ -Aᵧ * zz * √I / (√I_one + √I) / √I_one + # NOTE: The last √I_one here is not in the paper but is needed to make the units balance.
                 (zz / (ParentScope(cation.z) + ParentScope(anion.z))) *
                 (F_cat / ParentScope(cation.z) + F_an / ParentScope(anion.z))

        # Supplemental equations after equations 7 and 8
        Y ~ ((ParentScope(cation.z) + ParentScope(anion.z)) / 2)^2 * ParentScope(anion.m) /
            I
        X ~ ((ParentScope(cation.z) + ParentScope(anion.z)) / 2)^2 * ParentScope(cation.m) /
            I
        # Equation 9
        logγ⁰ ~ zz * logΓ⁰
        # Equation 10
        logΓ⁰ ~ log(1 + B * (1 + 0.1I / I_one)^q - B) + logΓ⁺
        # Equation 11
        B ~ 0.75 - 0.065q
        # Equation 12
        logΓ⁺ ~ -0.5107√I / (√I_one + C * √I)
        # Equation 13
        C ~ 1 + 0.055q * exp(-0.023I^3 / I_one^3)

        # Equation 14
        if is_CaSO4
            logγ ~ log(1.0) # From Table 4 footnote a, CaSO4 has an activity coefficient of zero,
        # but that would mean that either the Ca or the SO4 concentration would have to be
        # zero, so we set it to 1 here.
        elseif is_KHSO4 || is_NH4HSO4 || is_NaHSO4  # From Table 4 footnote b, c & d
            logγ ~ 0.5 * (ParentScope(salt1.logγ) + ParentScope(salt2.logγ) -
                ParentScope(salt3.logγ))
        elseif is_NH43HSO42
            logγ ~ (ParentScope(salt1.logγ) * 3 + ParentScope(salt2.logγ)) * (1 / 5) # From Table 4 footnote e
        else
            logγ ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀ -
                   (0.125 - c_1 * (T - T₀₂)) * A
        end

        # Equation in text below Equation 14
        A ~ -((0.41√I / (√I_one + √I)) + 0.039(I / I_one)^0.92)

        # Activity (Section 2.2)
        loga ~ ν_cation * log(ParentScope(cation.m) / m_one) +
               ν_anion * log(ParentScope(anion.m) / m_one) + (ν_cation + ν_anion) * logγ
    end
end

@mtkmodel BinaryMolality begin
    @description """Water content for a binary aqueous aerosol solution using data from
    Fountoukis and Nenes (2007) Table 7, for use in Equation 16."""
    @constants begin
        k_0, [description = "polynomial fit coefficient"]
        k_1, [description = "polynomial fit coefficient"]
        k_2, [description = "polynomial fit coefficient"]
        k_3, [description = "polynomial fit coefficient"]
        k_4, [description = "polynomial fit coefficient"]
        k_5, [description = "polynomial fit coefficient"]
        k_6, [description = "polynomial fit coefficient"]
        m_one = 1.0, [unit = u"mol/kg", description = "unit molality"]
    end
    @parameters begin
        RH, [description = "Relative humidity"]
    end
    @variables begin
        m_aw(t), [description = "molality of the binary solution", unit = u"mol/kg"]
    end
    @equations begin
        m_aw ~ (k_0 + k_1 * RH + k_2 * RH^2 + k_3 * RH^3 + k_4 * RH^4 + k_5 * RH^5) * m_one
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

        m_one = 1.0, [unit = u"mol/kg", description = "unit molality"]
    end
    @parameters begin
        T = 293.15, [description = "Temperature", unit = u"K"]
        RH = 0.3, [description = "Relative humidity (0-1)"]
    end
    @variables begin
        #! format: off
        I(t), [description = "Ionic strength", unit = u"mol/kg", guess=1.0, state_priority=10]
        W(t), [description = "Aerosol water content (per m^3 air)", unit = u"kg/m^3", guess=1e-6, state_priority=10]

        F_Ca(t), [description = "Activity contribution from Ca-containing salts"]
        F_K(t), [description = "Activity contribution from K-containing salts"]
        F_Mg(t), [description = "Activity contribution from Mg-containing salts"]
        F_Na(t), [description = "Activity contribution from Na-containing salts"]
        F_NH4(t), [description = "Activity contribution from NH4-containing salts"]
        F_Cl(t), [description = "Activity contribution from Cl-containing salts"]
        F_H(t), [description = "Activity contribution from H-containing salts"]
        F_NO3(t), [description = "Activity contribution from NO3-containing salts"]
        F_SO4(t), [description = "Activity contribution from SO4-containing salts"]
        F_HSO4(t), [description = "Activity contribution from HSO4-containing salts"]

        Aᵧ_term(t), [description = "Debye-Hückel term used in Equation 7 and 8"]
    end
    @structural_parameters begin
        q0 = 0.0 # Assume q=0 for salts in Table 4 with no data for q.
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
        HSO4 = Ion(z = abs(-1))
        OH = Ion(z = abs(-1))

        # Neutral species
        NH3 = Ion(z = 0)
        HNO3_aq = Ion(z = 0)
        HCl_aq = Ion(z = 0)

        # Salts
        CaNO32 = Salt(cation = Ca, anion = NO3, drh = 0.4906, l_t = 509.4, q = 0.93, T=T)
        CaCl2 = Salt(cation = Ca, anion = Cl, drh = 0.2830, l_t = 551.1, q = 2.4, T=T)
        K2SO4 = Salt(cation = K, anion = SO4, drh = 0.9751, l_t = 35.6, q = -0.25, T=T)
        KNO3 = Salt(cation = K, anion = NO3, drh = 0.9248, l_t = NaN, q = -2.33, T=T)
        KCl = Salt(cation = K, anion = Cl, drh = 0.8426, l_t = 158.9, q = 0.92, T=T)
        MgSO4 = Salt(cation = Mg, anion = SO4, drh = 0.8613, l_t = -714.5, q = 0.15, T=T)
        MgNO32 = Salt(cation = Mg, anion = NO3, drh = 0.5400, l_t = 230.2, q = 2.32, T=T)
        MgCl2 = Salt(cation = Mg, anion = Cl, drh = 0.3284, l_t = 42.23, q = 2.90, T=T)
        NaCl = Salt(cation = Na, anion = Cl, drh = 0.7528, l_t = 25.0, q = 2.23, T=T)
        Na2SO4 = Salt(cation = Na, anion = SO4, drh = 0.9300, l_t = 80.0, q = -0.19, T=T)
        NaNO3 = Salt(cation = Na, anion = NO3, drh = 0.7379, l_t = 304.0, q = -0.39, T=T)
        NH42SO4 = Salt(cation = NH4, anion = SO4, drh = 0.7997, l_t = 80.0, q = -0.25, T=T)
        NH4NO3 = Salt(cation = NH4, anion = NO3, drh = 0.6183, l_t = 852.0, q = -1.15, T=T)
        NH4Cl = Salt(cation = NH4, anion = Cl, drh = 0.7710, l_t = 239.0, q = 0.82, T=T)
        H2SO4 = Salt(cation = H, anion = SO4, drh = 0.000, l_t = NaN, q = -0.1, T=T)
        HHSO4 = Salt(cation = H, anion = HSO4, drh = 0.000, l_t = NaN, q = 8.00, T=T)
        HNO3 = Salt(cation = H, anion = NO3, drh = NaN, l_t = NaN, q = 2.60, T=T)
        HCl = Salt(cation = H, anion = Cl, drh = NaN, l_t = NaN, q = 6.00, T=T)
        CaSO4 = Salt(cation = Ca, anion = SO4, drh = 0.9700, l_t = NaN, q = q0, T=T,
            is_CaSO4=true)
        KHSO4 = Salt(cation = K, anion = HSO4, drh = 0.8600, l_t = NaN, q = q0, T=T,
            is_KHSO4=true, salt1=HHSO4, salt2=KCl, salt3=HCl)
        NH4HSO4 = Salt(cation = NH4, anion = HSO4, drh = 0.4000, l_t = 384.0, q = q0, T=T,
            is_NH4HSO4=true, salt1=HHSO4, salt2=NH4Cl, salt3=HCl)
        NaHSO4 = Salt(cation = Na, anion = HSO4, drh = 0.5200, l_t = -45.0, q = q0, T=T,
            is_NaHSO4=true, salt1=HHSO4, salt2=NaCl, salt3=HCl)
        NH43HSO42 = Salt(cation = NH4, anion = HSO4, drh = 0.6900, l_t = 186.0, q = q0, T=T,
            is_NH43HSO42=true, salt1=NH42SO4, salt2=NH4HSO4)

        # Water content
        #! format: off
        maw_CaNO32 = BinaryMolality(k_0=36.356, k_1=-165.66, k_2=447.46, k_3=-673.55, k_4=510.91, k_5=-155.56, k_6=0)
        maw_CaCl2 = BinaryMolality(k_0=20.847, k_1=-97.599, k_2=273.220, k_3=-422.120, k_4=331.160, k_5=-105.450, k_6=0)
        maw_KHSO4 = BinaryMolality(k_0=1.061, k_1=-0.101, k_2=1.579e-2, k_3=-1.950e-3, k_4=9.515e-5, k_5=-1.547e-6, k_6=0)
        maw_K2SO4 = BinaryMolality(k_0=1061.51, k_1=-4748.97, k_2=8096.16, k_3=-6166.16, k_4=1757.47, k_5=0, k_6=0)
        maw_KNO3 = BinaryMolality(k_0=1.2141e4, k_1=-5.1173e4, k_2=8.1252e4, k_3=-5.7527e4, k_4=1.5305e4, k_5=0, k_6=0)
        maw_KCl = BinaryMolality(k_0=179.721, k_1=-721.266, k_2=1161.03, k_3=-841.479, k_4=221 / 943, k_5=0, k_6=0)
        maw_MgSO4 = BinaryMolality(k_0=-0.778, k_1=177.74, k_2=-719.79, k_3=1174.6, k_4=-863.44, k_5=232.31, k_6=0)
        maw_MgNO32 = BinaryMolality(k_0=12.166, k_1=-16.154, k_2=0, k_3=10.886, k_4=0, k_5=-6.815, k_6=0)
        maw_MgCl2 = BinaryMolality(k_0=11.505, k_1=-26.518, k_2=34.937, k_3=-19.829, k_4=0, k_5=0, k_6=0)
        maw_NaNO3 = BinaryMolality(k_0=0.9988, k_1=-2.6947e-2, k_2=1.9610e-4, k_3=2.8154e-5, k_4=6.1359e-7, k_5=0, k_6=0)
        maw_NaHSO4 = BinaryMolality(k_0=1.0614, k_1=-0.1014, k_2=1.5796e-2, k_3=-1.9501e-3, k_4=9.5147e-5, k_5=-1.5473e-6, k_6=0)
        maw_NaCl = BinaryMolality(k_0=1.0084, k_1=-4.9390e-2, k_2=8.888e-3, k_3=-2.1570e-3, k_4=1.6170e-4, k_5=1.99e-6, k_6=-1.142e-7)
        maw_Na2SO4 = BinaryMolality(k_0=1.0052, k_1=-6.4840e-2, k_2=3.519e-2, k_3=-1.3190e-2, k_4=1.9250e-3, k_5=-1.224e-4, k_6=2.87e-6)
        maw_NH42SO4 = BinaryMolality(k_0=0.9968, k_1=-2.9690e-2, k_2=1.735e-5, k_3=-3.2530e-4, k_4=3.5710e-5, k_5=-9.7870e-7, k_6=0)
        maw_NH4Cl = BinaryMolality(k_0=0.9968, k_1=-2.6110e-2, k_2=-1.5990e-3, k_3=1.3550e-4, k_4=-2.317e-6, k_5=-1.113e-8, k_6=0)
        maw_NH4NO3 = BinaryMolality(k_0=1.0053, k_1=-2.4991e-2, k_2=4.4688e-4, k_3=1.6453e-5, k_4=-3.8940e-7, k_5=-4.7668e-8, k_6=1.3753e-9)
        maw_NH4HSO4 = BinaryMolality(k_0=1.0261, k_1=-4.9766e-2, k_2=3.2757e-3, k_3=-2.4477e-4, k_4=1.0766e-5, k_5=-1.8329e-7, k_6=0)
        maw_NH43HSO42 = BinaryMolality(k_0=1.0088, k_1=-5.3730e-2, k_2=1.4201e-3, k_3=-9.2484e-4, k_4=2.2796e-4, k_5=-1.5445e-5, k_6=0)
        #! format: on
    end
    @equations begin
        Aᵧ_term ~ Aᵧ * √I / (√I_one + √I) / √I_one # NOTE: The last √I_one here is not in the paper but is needed to make the units balance.
        # Equation 7
        F_NH4 ~ sum([s.Y * s.logγ⁰ for s in [NH4NO3, NH4Cl, NH4HSO4, NH42SO4, NH43HSO42]]) +
                Aᵧ_term *
                sum([s.zz * s.Y for s in [NH4NO3, NH4Cl, NH4HSO4, NH42SO4, NH43HSO42]])
        F_Na ~ sum([s.Y * s.logγ⁰ for s in [NaCl, Na2SO4, NaNO3, NaHSO4]]) +
               Aᵧ_term * sum([s.zz * s.Y for s in [NaCl, Na2SO4, NaNO3, NaHSO4]])
        F_H ~ sum([s.Y * s.logγ⁰ for s in [H2SO4, HCl, HNO3]]) +
              Aᵧ_term * sum([s.zz * s.Y for s in [H2SO4, HCl, HNO3]])
        F_Ca ~ sum([s.Y * s.logγ⁰ for s in [CaNO32, CaCl2, CaSO4]]) +
               Aᵧ_term * sum([s.zz * s.Y for s in [CaNO32, CaCl2, CaSO4]])
        F_K ~ sum([s.Y * s.logγ⁰ for s in [KHSO4, K2SO4, KNO3, KCl]]) +
              Aᵧ_term * sum([s.zz * s.Y for s in [KHSO4, K2SO4, KNO3, KCl]])
        F_Mg ~ sum([s.Y * s.logγ⁰ for s in [MgSO4, MgNO32, MgCl2]]) +
               Aᵧ_term * sum([s.zz * s.Y for s in [MgSO4, MgNO32, MgCl2]])

        # Equation 8
        F_Cl ~ sum([s.X * s.logγ⁰ for s in [NaCl, KCl, MgCl2, CaCl2, NH4Cl, HCl]]) +
               Aᵧ_term *
               sum([s.zz * s.X for s in [NaCl, KCl, MgCl2, CaCl2, NH4Cl, HCl]])
        F_NO3 ~ sum([s.X * s.logγ⁰ for s in [NaNO3, KNO3, MgNO32, CaNO32, NH4NO3, HNO3]]) +
                Aᵧ_term *
                sum([s.zz * s.X for s in [NaNO3, KNO3, MgNO32, CaNO32, NH4NO3, HNO3]])
        F_SO4 ~ sum([s.X * s.logγ⁰ for s in [Na2SO4, K2SO4, MgSO4, CaSO4, NH42SO4, H2SO4]]) +
                Aᵧ_term *
                sum([s.zz * s.X for s in [Na2SO4, K2SO4, MgSO4, CaSO4, NH42SO4, H2SO4]])
        F_HSO4 ~ sum([s.X * s.logγ⁰ for s in [KHSO4, HHSO4, NaHSO4, NH4HSO4, NH43HSO42]]) +
                 Aᵧ_term *
                 sum([s.zz * s.X for s in [KHSO4, HHSO4, NaHSO4, NH4HSO4, NH43HSO42]])

        CaNO32.F_cat ~ F_Ca
        CaNO32.F_an ~ F_NO3
        CaCl2.F_cat ~ F_Ca
        CaCl2.F_an ~ F_Cl
        CaSO4.F_cat ~ F_Ca
        CaSO4.F_an ~ F_SO4
        KHSO4.F_cat ~ F_K
        KHSO4.F_an ~ F_HSO4
        K2SO4.F_cat ~ F_K
        K2SO4.F_an ~ F_SO4
        KNO3.F_cat ~ F_K
        KNO3.F_an ~ F_NO3
        KCl.F_cat ~ F_K
        KCl.F_an ~ F_Cl
        MgSO4.F_cat ~ F_Mg
        MgSO4.F_an ~ F_SO4
        MgNO32.F_cat ~ F_Mg
        MgNO32.F_an ~ F_NO3
        MgCl2.F_cat ~ F_Mg
        MgCl2.F_an ~ F_Cl
        NaCl.F_cat ~ F_Na
        NaCl.F_an ~ F_Cl
        Na2SO4.F_cat ~ F_Na
        Na2SO4.F_an ~ F_SO4
        NaNO3.F_cat ~ F_Na
        NaNO3.F_an ~ F_NO3
        NH42SO4.F_cat ~ F_NH4
        NH42SO4.F_an ~ F_SO4
        NH4NO3.F_cat ~ F_NH4
        NH4NO3.F_an ~ F_NO3
        NH4Cl.F_cat ~ F_NH4
        NH4Cl.F_an ~ F_Cl
        NH4HSO4.F_cat ~ F_NH4
        NH4HSO4.F_an ~ F_HSO4
        NaHSO4.F_cat ~ F_Na
        NaHSO4.F_an ~ F_HSO4
        NH43HSO42.F_cat ~ F_NH4
        NH43HSO42.F_an ~ F_HSO4
        H2SO4.F_cat ~ F_H
        H2SO4.F_an ~ F_SO4
        HHSO4.F_cat ~ F_H
        HHSO4.F_an ~ F_HSO4
        HNO3.F_cat ~ F_H
        HNO3.F_an ~ F_NO3
        HCl.F_cat ~ F_H
        HCl.F_an ~ F_Cl


        # Ionic strength (equation in text below Equation 8)
        # We take the absolute value to avoid numerical issues as negative concentrations
        # sometimes occur during DAE initialization.
        I ~ abs(0.5 * sum([ion.m * ion.z^2
                     for ion in [
            NH4, Na, H, Ca, K, Mg, Cl, NO3, SO4, HSO4, OH]]))
        CaNO32.I ~ I
        CaCl2.I ~ I
        CaSO4.I ~ I
        KHSO4.I ~ I
        K2SO4.I ~ I
        KNO3.I ~ I
        KCl.I ~ I
        MgSO4.I ~ I
        MgNO32.I ~ I
        MgCl2.I ~ I
        NaCl.I ~ I
        Na2SO4.I ~ I
        NaNO3.I ~ I
        NH42SO4.I ~ I
        NH4NO3.I ~ I
        NH4Cl.I ~ I
        NH4HSO4.I ~ I
        NaHSO4.I ~ I
        NH43HSO42.I ~ I
        H2SO4.I ~ I
        HHSO4.I ~ I
        HNO3.I ~ I
        HCl.I ~ I


        # Water content (Section 2.3)
        W ~ CaNO32.M / maw_CaNO32.m_aw +
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
            NH43HSO42.M / maw_NH43HSO42.m_aw

        maw_CaNO32.RH ~ RH
        maw_CaCl2.RH ~ RH
        maw_KHSO4.RH ~ RH
        maw_K2SO4.RH ~ RH
        maw_KNO3.RH ~ RH
        maw_KCl.RH ~ RH
        maw_MgSO4.RH ~ RH
        maw_MgNO32.RH ~ RH
        maw_MgCl2.RH ~ RH
        maw_NaNO3.RH ~ RH
        maw_NaHSO4.RH ~ RH
        maw_NaCl.RH ~ RH
        maw_Na2SO4.RH ~ RH
        maw_NH42SO4.RH ~ RH
        maw_NH4Cl.RH ~ RH
        maw_NH4NO3.RH ~ RH
        maw_NH4HSO4.RH ~ RH
        maw_NH43HSO42.RH ~ RH

        H2SO4.M ~ 0.0 # The first dissociation of H2SO4 is assumed to be complete (Section 3.3)
        CaSO4.M ~ 0.0 # From Table 4 footnote a, all of CaSO4 precipitates to the solid phase.

        # FIXME(CT): I can't figure out how to get the aqueous salt mass balance
        # to work consistently, so here is a had to just evenly split the mass
        # of each cation between its corresponding salts.
        # See below for other aborted attempts.
        NH4NO3.M ~ NH4.m * W / 5
        NH4Cl.M ~ NH4.m * W / 5
        NH4HSO4.M ~ NH4.m * W / 5
        2NH42SO4.M ~ NH4.m * W / 5
        3NH43HSO42.M ~ NH4.m * W / 5

        NaCl.M ~ Na.m * W / 4
        2Na2SO4.M ~ Na.m * W / 4
        NaNO3.M ~ Na.m * W / 4
        NaHSO4.M ~ Na.m * W / 4

        HCl.M ~ H.m * W / 3
        HNO3.M ~ H.m * W / 3
        HHSO4.M ~ H.m * W / 3

        CaNO32.M ~ Ca.m * W / 2
        CaCl2.M ~ Ca.m * W / 2

        KHSO4.M ~ K.m * W / 4
        2K2SO4.M ~ K.m * W / 4
        KNO3.M ~ K.m * W / 4
        KCl.M ~ K.m * W / 4

        MgSO4.M ~ Mg.m * W / 3
        MgNO32.M ~ Mg.m * W / 3
        #MgCl2.M ~ Mg.m * W / 3 # Leave one degree of freedom to set overall molarity level.
        # FIXME(CT): This can result in negative MgCl2.M, which can effect water concentration.

        # Mass balance
        # NH4.m * W ~ sum([NH4NO3.M, NH4Cl.M, NH4HSO4.M, 2NH42SO4.M, 3NH43HSO42.M])
        # Na.m * W ~ sum([NaCl.M, 2Na2SO4.M, NaNO3.M, NaHSO4.M])
        # H.m * W ~ sum([2H2SO4.M, HCl.M, HNO3.M, HHSO4.M])
        # Ca.m * W ~ sum([CaNO32.M, CaCl2.M, CaSO4.M])
        # K.m * W ~ sum([KHSO4.M, 2K2SO4.M, KNO3.M, KCl.M])
        # Mg.m * W ~ sum([MgSO4.M, MgNO32.M, MgCl2.M])
        # Cl.m * W ~ sum([NaCl.M, KCl.M, 2MgCl2.M, 2CaCl2.M, NH4Cl.M, HCl.M])
        # NO3.m * W ~ sum([NaNO3.M, KNO3.M, 2MgNO32.M, 2CaNO32.M, NH4NO3.M, HNO3.M])
        # SO4.m * W ~ sum([Na2SO4.M, K2SO4.M, MgSO4.M, CaSO4.M, NH42SO4.M, H2SO4.M])
        # HSO4.m * W ~ sum([KHSO4.M, NaHSO4.M, NH4HSO4.M, 2NH43HSO42.M, HHSO4.M])

        # NH43HSO42.M / HHSO4.M ~ (a_NH43HSO42 / a_HHSO4) / m_one^3
        # CaCl2.M / NaCl.M ~ (a_CaCl2 / a_NaCl) / m_one
        # K2SO4.M / KCl.M ~ (a_K2SO4 / a_KCl) / m_one
        # MgNO32.M / KNO3.M ~ (a_MgNO32 / a_KNO3) / m_one
        # Na2SO4.M / NaHSO4.M ~ (a_Na2SO4 / a_NaHSO4) / m_one
        # NH42SO4.M / NH4NO3.M ~ (a_NH42SO4 / a_NH4NO3) / m_one
        # MgCl2.M / NH4Cl.M ~ (a_MgCl2 / a_NH4Cl) / m_one
        # MgSO4.M / KHSO4.M ~ (a_MgSO4 / a_KHSO4)
        # NH4HSO4.M / HCl.M ~ (a_NH4HSO4 / a_HCl)
        # NaNO3.M / HNO3.M ~ (a_NaNO3 / a_HNO3)
        # Leftover: CaNO32

        # 0 ~ min(CaNO32.M, CaCl2.M, CaSO4.M, KHSO4.M, K2SO4.M, KNO3.M, KCl.M, MgSO4.M,
        # MgNO32.M, MgCl2.M, NaCl.M, Na2SO4.M, NaNO3.M, NH42SO4.M, NH4NO3.M,
        # NH4Cl.M, NH4HSO4.M, NaHSO4.M, NH43HSO42.M, H2SO4.M, HHSO4.M, HNO3.M, HCl.M)

        # Second mass balance to promote non-negativity
        #         sum([NH4.m, Na.m, H.m, Ca.m, K.m, Mg.m,
        #         Cl.m, NO3.m, SO4.m, HSO4.m, OH.m, NH3.m, HNO3_aq.m, HCl_aq.m, H.m, OH.m
        # ]) ~ sum(abs.([NH4.m, Na.m, H.m, Ca.m, K.m, Mg.m,
        #             Cl.m, NO3.m, SO4.m, HSO4.m, OH.m, NH3.m, HNO3_aq.m, HCl_aq.m, H.m, OH.m]))

        #         sum([NH4NO3.M, NH4Cl.M, NH4HSO4.M, 2NH42SO4.M, 3NH43HSO42.M,
        #         NaCl.M, 2Na2SO4.M, NaNO3.M, NaHSO4.M,
        #         2H2SO4.M, HCl.M, HNO3.M, KHSO4.M, HHSO4.M,
        #         CaNO32.M, CaCl2.M, CaSO4.M,
        #         KHSO4.M, 2K2SO4.M, KNO3.M, KCl.M,
        #         MgSO4.M, MgNO32.M, MgCl2.M
        # ]) ~ sum(abs.([NH4NO3.M, NH4Cl.M, NH4HSO4.M, 2NH42SO4.M, 3NH43HSO42.M,
        #             NaCl.M, 2Na2SO4.M, NaNO3.M, NaHSO4.M,
        #             2H2SO4.M, HCl.M, HNO3.M, KHSO4.M, HHSO4.M,
        #             CaNO32.M, CaCl2.M, CaSO4.M,
        #             KHSO4.M, 2K2SO4.M, KNO3.M, KCl.M,
        #             MgSO4.M, MgNO32.M, MgCl2.M]))

        # Charge balance
        0 ~ sum([i.m * i.z for i in [NH4, Na, H, Ca, K, Mg]]) -
            sum([i.m * i.z for i in [Cl, NO3, SO4, HSO4, OH]])
    end
end
