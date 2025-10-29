@mtkmodel Ion begin
    @description "An aqueous ion."
    @structural_parameters begin
        m_eq_guess = 1 # Guess for equilibrium molality
    end
    @parameters begin
        z, [description = "Valence (charge) of the ion"]
    end
    @variables begin
        #! format: off
        m_eq(t), [description = "Molality of ion in water at equilibrium", unit = u"mol/kg", guess = m_eq_guess, state_priority=10]
        m_aq(t), [description = "Molality of ion in water", unit = u"mol/kg", guess = 1]
        M(t), [description = "Molarity of ion in air", unit = u"mol/m^3", guess = 1e-8]
        W_eq(t), [description = "Aerosol water content in air at equilibrium", unit = u"kg/m^3", guess = 1e-6]
        W(t), [description = "Aerosol water content in air", unit = u"kg/m^3", guess = 1e-6]
        #! format: on
    end
    @equations begin
        M ~ m_eq * abs(W_eq)
        m_aq ~ M / abs(W)
    end
end

@mtkmodel Salt begin
    @description """
An aqueous salt comprised of a cation, an anion, and an activity parameter (q).
q values are given in Table 4 of Fountoukis and Nenes (2007).
"""
    @structural_parameters begin
        is_KHSO4 = false # Is this salt KHSO4?
        is_NH4HSO4 = false # Is this salt NH4HSO4?
        is_NaHSO4 = false # Is this salt NaHSO4?
        is_NH43HSO42 = false # Is this salt NH43HSO42?
        salt1 = nothing
        salt2 = nothing
        salt3 = nothing
        m_eq_guess = 1 # Guess for equilibrium molality
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
        z_cation = 1, [description = "Absolute value of the cation charge"]
        z_anion = 1, [description = "Absolute value of the anion charge"]

        # Derived constants
        B
        zz, [description = "Product of the absolute values of the cation and anion charges"]

        # Unit conversions
        I_one = 1, [unit = u"mol/kg", description = "An ionic strength of 1"]
        m_one = 1, [unit = u"mol/kg", description = "A molality of 1"]
    end
    @parameters begin
        T = 293.15, [description = "Temperature", unit = u"K"]
        RH = 0.3, [description = "Relative humidity (0-1)"]
    end
    @variables begin
        #! format: off
        m_eq(t), [description = "molality of the salt in water at equilibrium", unit = u"mol/kg", guess = m_eq_guess, state_priority=10]
        m_aq(t), [description = "molality of the salt in water", unit = u"mol/kg", guess = 1]
        M_eq(t), [description = "Molarity of the salt in air at equilibrium", unit = u"mol/m^3", guess = 1e-8]
        M_aq(t), [description = "Molarity of the salt in air", unit = u"mol/m^3", guess = 1e-8]
        M_precip(t), [description = "Precipitated solid", unit = u"mol/m^3", guess=1e-10]
        M_dissolved(t), [description = "Dissolved precipitate", unit = u"mol/m^3", guess=1e-10]
        M_salt(t), [description = "Salt concentration (either precipitated or dissolved)", unit = u"mol/m^3", guess=1e-10]
        deliquesced(t), [description = "Whether the salt is deliquesced (1) or not (0)"]
        W_eq(t), [description = "Aerosol water content in air associated with equilibrium concentration", unit = u"kg/m^3", guess = 1e-6]
        W(t), [description = "Aerosol water content in air", unit = u"kg/m^3", guess = 1e-6]
        X(t)
        Y(t)
        I(t), [description = "Ionic strength", unit = u"mol/kg", state_priority=-10]
        logγ⁰(t), [description = "Log of the mean ionic activity coefficient for the single solute solution"]
        logγₜ₀(t), [description = "Log of the multi-component activity coefficient at 298.15K", state_priority=-10]
        logγ(t), [description = "Log of the activity coefficient"]
        #! format: on
        logΓ⁰(t)
        logΓ⁺(t)
        C(t)
        F_cat(t), [description = "Activity contribution from the cation"]
        F_an(t), [description = "Activity contribution from the anion"]
        A(t), [description = "Parameter used in Equation 14"]
        loga_eq(t), [description = "Log of the activity of the salt at equilibrium"]
        loga_aq(t), [description = "Log of the aqueous activity of the salt"]
    end
    @equations begin
        M_eq ~ m_eq * abs(W_eq)
        m_aq ~ M_aq / abs(W)
        zz ~ z_cation * z_anion

        # Equation 6
        logγₜ₀ ~ -Aᵧ * zz * √I / (√I_one + √I) / √I_one + # NOTE: The last √I_one here is not in the paper but is needed to make the units balance.
                 (zz / (z_cation + z_anion)) *
                 (F_cat / z_cation + F_an / z_anion)

        # Supplemental equations after equations 7 and 8
        Y ~ ((z_cation + z_anion) / 2)^2 * ν_anion * m_eq / I
        X ~ ((z_cation + z_anion) / 2)^2 * ν_cation * m_eq / I
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
        if is_KHSO4 || is_NH4HSO4 || is_NaHSO4  # From Table 4 footnote b, c & d
            logγ ~ 0.5 * (ParentScope(salt1.logγ) + ParentScope(salt2.logγ) -
                    ParentScope(salt3.logγ))
        elseif is_NH43HSO42
            logγ ~ (ParentScope(salt1.logγ) * 3 + ParentScope(salt2.logγ)) * (1 / 5) # From Table 4 footnote e
        else
            logγ ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀ - (0.125 - c_1 * (T - T₀₂)) * A
        end

        # Equation in text below Equation 14
        A ~ -((0.41√I / (√I_one + √I)) + 0.039(I / I_one)^0.92)

        # Activity (Section 2.2)
        loga_eq ~ ν_cation * log(m_eq / m_one) + ν_anion * log(m_eq / m_one) +
                  (ν_cation + ν_anion) * logγ
        loga_aq ~ ν_cation * log(m_aq / m_one) + ν_anion * log(m_aq / m_one) +
                  (ν_cation + ν_anion) * logγ

        deliquesced ~  (tanh((RH - drh) * 30) + 1) / 2
        M_precip ~ M_salt * (1 - deliquesced)
         M_dissolved ~ M_salt * deliquesced
        M_eq ~ M_aq - M_dissolved
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
        W_eq(t), [description = "Aerosol water content in air associated with equilibrium", unit = u"kg/m^3", guess=1e-6, state_priority=10]
        W(t), [description = "Aerosol water content in air", unit = u"kg/m^3", guess=1e-6]

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
        F_OH(t), [description = "Activity contribution from OH-containing salts"]

        Aᵧ_term(t), [description = "Debye-Hückel term used in Equation 7 and 8"]
    end
    @structural_parameters begin
        q0 = 0.0 # Assume q=0 for salts in Table 4 with no data for q.
    end
    @components begin
        # Cations
        Na = Ion(z = 1, m_eq_guess=160.0)
        H = Ion(z = 1, m_eq_guess=460.0)
        Ca = Ion(z = 2, m_eq_guess=12.0)
        K = Ion(z = 1, m_eq_guess=5.9)
        Mg = Ion(z = 2, m_eq_guess=130.0)

        # Anions
        OH = Ion(z = abs(-1), m_eq_guess=13.0)

        # Neutral species
        NH3 = Ion(z = 0, m_eq_guess=5.8e-5)
        HNO3_aq = Ion(z = 0, m_eq_guess=0.056)
        HCl_aq = Ion(z = 0, m_eq_guess=0.0018)

        # Salts
        CaNO32 = Salt(z_cation = 2, ν_anion=2, drh = 0.4906, l_t = 509.4, q = 0.93, T=T,
             RH=RH,m_eq_guess=5.6)
        CaCl2 = Salt(z_cation = 2, ν_anion=2, drh = 0.2830, l_t = 551.1, q = 2.4, T=T,
            RH=RH,m_eq_guess=6.3)
        CaSO4 = Salt(z_cation = 2, z_anion = 2, drh = 0.9700, l_t = NaN, q = q0, T=T,
             RH=RH,m_eq_guess=0.065)
        K2SO4 = Salt(ν_cation = 2, z_anion = 2, drh = 0.9751, l_t = 35.6, q = -0.25, T=T,
             RH=RH,m_eq_guess=0.42)
        KNO3 = Salt(z_cation = 1, z_anion = 1, drh = 0.9248, l_t = NaN, q = -2.33, T=T,
             RH=RH,m_eq_guess=0.11)
        KCl = Salt(z_cation = 1, z_anion = 1, drh = 0.8426, l_t = 158.9, q = 0.92, T=T,
             RH=RH,m_eq_guess=0.011)
        MgSO4 = Salt(z_cation = 2, z_anion = 2, drh = 0.8613, l_t = -714.5, q = 0.15, T=T,
            RH=RH,m_eq_guess=0.42)
        MgNO32 = Salt(z_cation = 2, ν_anion = 2, drh = 0.5400, l_t = 230.2, q = 2.32, T=T,
            RH=RH,m_eq_guess=52.0)
        MgCl2 = Salt(z_cation = 2, ν_anion = 2, drh = 0.3284, l_t = 42.23, q = 2.90, T=T,
            RH=RH,m_eq_guess=82.0)
        NaCl = Salt(z_cation = 1, z_anion = 1, drh = 0.7528, l_t = 25.0, q = 2.23, T=T,
            RH=RH,m_eq_guess=0.022)
        Na2SO4 = Salt(ν_cation = 2, z_anion = 2, drh = 0.9300, l_t = 80.0, q = -0.19, T=T,
            RH=RH,m_eq_guess=1.3)
        NaNO3 = Salt(z_cation = 1, z_anion = 1, drh = 0.7379, l_t = 304.0, q = -0.39, T=T,
             RH=RH,m_eq_guess=0.39)
        NH42SO4 = Salt(ν_cation = 2, z_anion = 2, drh = 0.7997, l_t = 80.0, q = -0.25, T=T,
             RH=RH,m_eq_guess=2.1)
        NH4NO3 = Salt(z_cation = 1, z_anion = 1, drh = 0.6183, l_t = 852.0, q = -1.15, T=T,
         RH=RH)
        NH4Cl = Salt(z_cation = 1, z_anion = 1, drh = 0.7710, l_t = 239.0, q = 0.82, T=T,
         RH=RH)
        HHSO4 = Salt(z_cation = 1, z_anion = 1, drh = 0.000, l_t = NaN, q = 8.00, T=T,
         RH=RH)
        H2SO4 = Salt(z_cation = 1, z_anion = 1, drh = 0.000, l_t = NaN, q = -0.1, T=T,
             RH=RH,m_eq_guess=160.0)
        HNO3 = Salt(z_cation = 1, z_anion = 1, drh = 0, l_t = NaN, q = 2.60, T=T,
             RH=RH,m_eq_guess=120.0)
        HCl = Salt(z_cation = 1, z_anion = 1, drh = 0, l_t = NaN, q = 6.00, T=T,
             RH=RH,m_eq_guess=180.0)
        KHSO4 = Salt(z_cation = 1, z_anion = 1, drh = 0.8600, l_t = NaN, q = q0, T=T,
             RH=RH,is_KHSO4=true, salt1=HHSO4, salt2=KCl, salt3=HCl, m_eq_guess=4.9)
        NH4HSO4 = Salt(z_cation = 1, z_anion = 1, drh = 0.4000, l_t = 384.0, q = q0, T=T,
             RH=RH,is_NH4HSO4=true, salt1=HHSO4, salt2=NH4Cl, salt3=HCl, m_eq_guess=1.2)
        NaHSO4 = Salt(z_cation = 1, z_anion = 1, drh = 0.5200, l_t = -45.0, q = q0, T=T,
             RH=RH,is_NaHSO4=true, salt1=HHSO4, salt2=NaCl, salt3=HCl, m_eq_guess=160.0)
        NH43HSO42 = Salt(ν_cation = 3, ν_anion = 2, drh = 0.6900, l_t = 186.0, q = q0, T=T,
             RH=RH,is_NH43HSO42=true, salt1=NH42SO4, salt2=NH4HSO4, m_eq_guess=2.7)

        # Species that are not in the paper. Assume q=0
        HSO4_dissociated = Salt(z_anion = abs(-2), drh = 0, l_t = NaN, q = 0.00, T=T,
            m_eq_guess=7.0) # H + SO4
        NH3_dissociated = Salt(drh = 0, l_t = NaN, q = 0.00, T=T,
            m_eq_guess=13.0) # NH4 + OH
        H2O_dissociated = Salt(drh = 0, l_t = NaN, q = 0.00, T=T,
            m_eq_guess=3.2e-7) # H + OH

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
        F_H ~ sum([s.Y * s.logγ⁰ for s in [HCl, HNO3]]) +
              Aᵧ_term * sum([s.zz * s.Y for s in [HCl, HNO3]])
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
        F_SO4 ~ sum([s.X * s.logγ⁰ for s in [Na2SO4, K2SO4, MgSO4, CaSO4, NH42SO4]]) +
                Aᵧ_term *
                sum([s.zz * s.X for s in [Na2SO4, K2SO4, MgSO4, CaSO4, NH42SO4]])
        F_HSO4 ~ sum([s.X * s.logγ⁰ for s in [KHSO4, NaHSO4, NH4HSO4, NH43HSO42]]) +
                 Aᵧ_term *
                 sum([s.zz * s.X for s in [KHSO4, NaHSO4, NH4HSO4, NH43HSO42]])
        F_OH ~ sum([s.X * s.logγ⁰ for s in [H2O_dissociated, NH3_dissociated]]) +
               Aᵧ_term *
               sum([s.zz * s.X for s in [H2O_dissociated, NH3_dissociated]])

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
        HHSO4.F_cat ~ F_H
        HHSO4.F_an ~ F_HSO4
        H2SO4.F_cat ~ F_H
        H2SO4.F_an ~ F_HSO4
        HNO3.F_cat ~ F_H
        HNO3.F_an ~ F_NO3
        HCl.F_cat ~ F_H
        HCl.F_an ~ F_Cl
        NH3_dissociated.F_cat ~ F_NH4
        NH3_dissociated.F_an ~ F_OH
        H2O_dissociated.F_cat ~ F_H
        H2O_dissociated.F_an ~ F_OH
        HSO4_dissociated.F_cat ~ F_H
        HSO4_dissociated.F_an ~ F_SO4

        # Ionic strength (equation in text below Equation 8)
        # We take the absolute value to avoid numerical issues as negative concentrations
        # sometimes occur during DAE initialization.
        # We also clip the ionic strength at 10 to avoid unrealistic values at low RH.
        I ~ abs(0.5 * (sum([salt.m_eq * salt.ν_cation * salt.z_cation^2 +
                      salt.m_eq * salt.ν_anion * salt.z_anion^2
                      for salt in [
                     CaNO32, CaCl2, CaSO4, KHSO4, K2SO4, KNO3, KCl, MgSO4, MgNO32,
                     MgCl2, NaCl, Na2SO4, NaNO3, NH42SO4, NH4NO3, NH4Cl, NH4HSO4,
                     NaHSO4, NH43HSO42]]) +
                 H.m_eq * H.z^2 + OH.m_eq * OH.z^2))

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
        HHSO4.I ~ I
        H2SO4.I ~ I
        HNO3.I ~ I
        HCl.I ~ I
        NH3_dissociated.I ~ I
        H2O_dissociated.I ~ I
        HSO4_dissociated.I ~ I

        # Water content (Section 2.3)
        W ~ CaNO32.M_aq / maw_CaNO32.m_aw +
            CaCl2.M_aq / maw_CaCl2.m_aw +
            KHSO4.M_aq / maw_KHSO4.m_aw +
            K2SO4.M_aq / maw_K2SO4.m_aw +
            KNO3.M_aq / maw_KNO3.m_aw +
            KCl.M_aq / maw_KCl.m_aw +
            MgSO4.M_aq / maw_MgSO4.m_aw +
            MgNO32.M_aq / maw_MgNO32.m_aw +
            MgCl2.M_aq / maw_MgCl2.m_aw +
            NaCl.M_aq / maw_NaCl.m_aw +
            Na2SO4.M_aq / maw_Na2SO4.m_aw +
            NaNO3.M_aq / maw_NaNO3.m_aw +
            NH42SO4.M_aq / maw_NH42SO4.m_aw +
            NH4NO3.M_aq / maw_NH4NO3.m_aw +
            NaHSO4.M_aq / maw_NaHSO4.m_aw +
            NH4Cl.M_aq / maw_NH4Cl.m_aw +
            NH4HSO4.M_aq / maw_NH4HSO4.m_aw +
            NH43HSO42.M_aq / maw_NH43HSO42.m_aw

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

        CaNO32.W_eq ~ W_eq
        CaCl2.W_eq ~ W_eq
        CaSO4.W_eq ~ W_eq
        KHSO4.W_eq ~ W_eq
        K2SO4.W_eq ~ W_eq
        KNO3.W_eq ~ W_eq
        KCl.W_eq ~ W_eq
        MgSO4.W_eq ~ W_eq
        MgNO32.W_eq ~ W_eq
        MgCl2.W_eq ~ W_eq
        NaCl.W_eq ~ W_eq
        Na2SO4.W_eq ~ W_eq
        NaNO3.W_eq ~ W_eq
        NH42SO4.W_eq ~ W_eq
        NH4NO3.W_eq ~ W_eq
        NH4Cl.W_eq ~ W_eq
        NH4HSO4.W_eq ~ W_eq
        NaHSO4.W_eq ~ W_eq
        NH43HSO42.W_eq ~ W_eq
        HHSO4.W_eq ~ W_eq
        H2SO4.W_eq ~ W_eq
        HNO3.W_eq ~ W_eq
        HCl.W_eq ~ W_eq

        Na.W_eq ~ W_eq
        H.W_eq ~ W_eq
        Ca.W_eq ~ W_eq
        K.W_eq ~ W_eq
        Mg.W_eq ~ W_eq
        HSO4_dissociated.W_eq ~ W_eq
        OH.W_eq ~ W_eq
        NH3.W_eq ~ W_eq
        NH3_dissociated.W_eq ~ W_eq
        HNO3_aq.W_eq ~ W_eq
        HCl_aq.W_eq ~ W_eq
        H2O_dissociated.W_eq ~ W_eq

        CaNO32.W ~ W
        CaCl2.W ~ W
        CaSO4.W ~ W
        KHSO4.W ~ W
        K2SO4.W ~ W
        KNO3.W ~ W
        KCl.W ~ W
        MgSO4.W ~ W
        MgNO32.W ~ W
        MgCl2.W ~ W
        NaCl.W ~ W
        Na2SO4.W ~ W
        NaNO3.W ~ W
        NH42SO4.W ~ W
        NH4NO3.W ~ W
        NH4Cl.W ~ W
        NH4HSO4.W ~ W
        NaHSO4.W ~ W
        NH43HSO42.W ~ W
        HHSO4.W ~ W
        H2SO4.W ~ W
        HNO3.W ~ W
        HCl.W ~ W

        Na.W ~ W
        H.W ~ W
        Ca.W ~ W
        K.W ~ W
        Mg.W ~ W
        HSO4_dissociated.W ~ W
        OH.W ~ W
        NH3.W ~ W
        NH3_dissociated.W ~ W
        HNO3_aq.W ~ W
        HCl_aq.W ~ W
        H2O_dissociated.W ~ W

        # NH4Cl and NH4NO3 precipitate directly from gas, so we assume
        # the aqueous concentration is zero.
        NH4Cl.m_eq ~ 1e-20 * m_one
        NH4NO3.m_eq ~ 1e-20 * m_one

        # HHSO4 is only used to calculate activity coefficients of other salts.
        HHSO4.m_eq ~ 1e-20 * m_one

        # Mass balance
        NH3_dissociated.M_eq ~ sum([NH4NO3.M_eq, NH4Cl.M_eq, NH4HSO4.M_eq, 2NH42SO4.M_eq,
            3NH43HSO42.M_eq])
        Na.M ~ sum([NaCl.M_eq, 2Na2SO4.M_eq, NaNO3.M_eq, NaHSO4.M_eq])
        H.M ~ sum([
            HCl.M_eq, H2SO4.M_eq, HSO4_dissociated.M_eq, H2O_dissociated.M_eq, HNO3.M_eq])
        Ca.M ~ sum([CaNO32.M_eq, CaCl2.M_eq, CaSO4.M_eq])
        K.M ~ sum([KHSO4.M_eq, 2K2SO4.M_eq, KNO3.M_eq, KCl.M_eq])
        Mg.M ~ sum([MgSO4.M_eq, MgNO32.M_eq, MgCl2.M_eq])
        HCl.M_eq ~ sum([NaCl.M_eq, KCl.M_eq, 2MgCl2.M_eq, 2CaCl2.M_eq, NH4Cl.M_eq])
        HNO3.M_eq ~ sum([NaNO3.M_eq, KNO3.M_eq, 2MgNO32.M_eq, 2CaNO32.M_eq, NH4NO3.M_eq])
        HSO4_dissociated.M_eq ~ sum([Na2SO4.M_eq, K2SO4.M_eq, MgSO4.M_eq, CaSO4.M_eq,
            NH42SO4.M_eq, NH43HSO42.M_eq])
        H2SO4.M_eq ~ sum([KHSO4.M_eq, NaHSO4.M_eq, NH4HSO4.M_eq, NH43HSO42.M_eq])
        OH.M ~ NH3_dissociated.M_eq + H2O_dissociated.M_eq
    end
end
