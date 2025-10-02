using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities

@mtkmodel Ion begin
    @description "An aqueous ion."
    @parameters begin
        z, [description = "Valence (charge) of the ion"]
    end
    @variables begin
        m(t) = 0.0, [description = "Molal concentration of the ion", unit = u"mol/kg"]
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
        cation = Ion(; z)
        anion = Ion(; z)
    end
    @constants begin
        # Salt properties
        drh, [description = "Deliquescence relative humidity at 298.15K"]
        l_t, [description = "Enthalpy term (-18/1000R L_s m_s)"]
        q, [description = "Kusik-Meissner Binary activity parameter"]

        # Derived constants
        B
        zz, [description = "Product of the absolute values of the cation and anion charges"]

        # Unit conversions
        I_one = 1, [unit = u"mol/kg", description = "An ionic strength of 1"]
    end
    @variables begin
        m(t), [description = "Molal concentration of the salt", unit = u"mol/kg"]
        X(t)
        Y(t)
        I(t), [description = "Ionic strength", unit = u"mol/kg"]
        logγ⁰(t), [description = "Log of the standard state activity coefficient"]
        Γ⁰(t), [unit = u"mol/kg"]
        Γ⁺(t)
        C(t)
    end
    @equations begin
        zz ~ ParentScope(cation.z) * ParentScope(anion.z)
        # Supplemental equations after equations 7 and 8
        Y ~ ((ParentScope(cation.z) + ParentScope(anion.z)) / 2)^2 * ParentScope(anion.m) /
            I
        X ~ ((ParentScope(cation.z) + ParentScope(anion.z)) / 2)^2 * ParentScope(cation.m) /
            I
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
    @variables begin
        RH(t), [description = "Relative humidity"]
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

        T₀₂ = 273.15, [unit = u"K", description = "Standard temperature 2"]
        c_1 = 0.005,
        [
            unit = u"K^-1",
            description = "Constant for Fountoukis and Nenes (2007) Eq. 14"
        ]
    end
    @parameters begin
        T = 298.15, [description = "Temperature", unit = u"K"]
        RH = 0.5, [description = "Relative humidity (0-1)"]
    end
    @variables begin
        I(t), [description = "Ionic strength", unit = u"mol/kg"]
        W(t), [description = "Aerosol water content (per m^3 air)", unit = u"kg/m^3"]

        #! format: off
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
        A(t), [description = "Parameter used in Equation 14"]

        logγₜ₀_CaNO32(t), [description = "Log of the standard state activity coefficient of CaNO3_2 at 298.15K"]
        logγₜ₀_CaCl2(t), [description = "Log of the standard state activity coefficient of CaCl2 at 298.15K"]
        logγₜ₀_K2SO4(t), [description = "Log of the standard state activity coefficient of K2SO4 at 298.15K"]
        logγₜ₀_KNO3(t), [description = "Log of the standard state activity coefficient of KNO3 at 298.15K"]
        logγₜ₀_KCl(t), [description = "Log of the standard state activity coefficient of KCl at 298.15K"]
        logγₜ₀_MgSO4(t), [description = "Log of the standard state activity coefficient of MgSO4 at 298.15K"]
        logγₜ₀_MgNO32(t), [description = "Log of the standard state activity coefficient of MgNO3_2 at 298.15K"]
        logγₜ₀_MgCl2(t), [description = "Log of the standard state activity coefficient of MgCl2 at 298.15K"]
        logγₜ₀_NaCl(t), [description = "Log of the standard state activity coefficient of NaCl at 298.15K"]
        logγₜ₀_Na2SO4(t), [description = "Log of the standard state activity coefficient of Na2SO4 at 298.15K"]
        logγₜ₀_NaNO3(t), [description = "Log of the standard state activity coefficient of NaNO3 at 298.15K"]
        logγₜ₀_NH42SO4(t), [description = "Log of the standard state activity coefficient of (NH4)2SO4 at 298.15K"]
        logγₜ₀_NH4NO3(t), [description = "Log of the standard state activity coefficient of NH4NO3 at 298.15K"]
        logγₜ₀_NH4Cl(t), [description = "Log of the standard state activity coefficient of NH4Cl at 298.15K"]
        logγₜ₀_H2SO4(t), [description = "Log of the standard state activity coefficient of H2SO4 at 298.15K"]
        logγₜ₀_HHSO4(t), [description = "Log of the standard state activity coefficient of HHSO4 at 298.15K"]
        logγₜ₀_HNO3_aqs(t), [description = "Log of the standard state activity coefficient of HNO3 (aqueous) at 298.15K"]
        logγₜ₀_HCl_aqs(t), [description = "Log of the standard state activity coefficient of HCl (aqueous) at 298.15K"]
        #! format: on

        logγ_CaNO32(t), [description = "Log of the activity coefficient of CaNO3_2"]
        logγ_CaCl2(t), [description = "Log of the activity coefficient of CaCl2"]
        logγ_K2SO4(t), [description = "Log of the activity coefficient of K2SO4"]
        logγ_KNO3(t), [description = "Log of the activity coefficient of KNO3"]
        logγ_KCl(t), [description = "Log of the activity coefficient of KCl"]
        logγ_MgSO4(t), [description = "Log of the activity coefficient of MgSO4"]
        logγ_MgNO32(t), [description = "Log of the activity coefficient of MgNO3_2"]
        logγ_MgCl2(t), [description = "Log of the activity coefficient of MgCl2"]
        logγ_NaCl(t), [description = "Log of the activity coefficient of NaCl"]
        logγ_Na2SO4(t), [description = "Log of the activity coefficient of Na2SO4"]
        logγ_NaNO3(t), [description = "Log of the activity coefficient of NaNO3"]
        logγ_NH42SO4(t), [description = "Log of the activity coefficient of (NH4)2SO4"]
        logγ_NH4NO3(t), [description = "Log of the activity coefficient of NH4NO3"]
        logγ_NH4Cl(t), [description = "Log of the activity coefficient of NH4Cl"]
        logγ_H2SO4(t), [description = "Log of the activity coefficient of H2SO4"]
        logγ_HHSO4(t), [description = "Log of the activity coefficient of HHSO4"]
        logγ_HNO3_aqs(t), [description = "Log of the activity coefficient of HNO3 (salt)"]
        logγ_HCl_aqs(t), [description = "Log of the activity coefficient of HCl (salt)"]

        γ_CaNO32(t), [description = "Activity coefficient of CaNO3_2"]
        γ_CaCl2(t), [description = "Activity coefficient of CaCl2"]
        γ_CaSO4(t), [description = "Activity coefficient of CaSO4"]
        γ_KHSO4(t), [description = "Activity coefficient of KHSO4"]
        γ_K2SO4(t), [description = "Activity coefficient of K2SO4"]
        γ_KNO3(t), [description = "Activity coefficient of KNO3"]
        γ_KCl(t), [description = "Activity coefficient of KCl"]
        γ_MgSO4(t), [description = "Activity coefficient of MgSO4"]
        γ_MgNO32(t), [description = "Activity coefficient of MgNO3_2"]
        γ_MgCl2(t), [description = "Activity coefficient of MgCl2"]
        γ_NaCl(t), [description = "Activity coefficient of NaCl"]
        γ_Na2SO4(t), [description = "Activity coefficient of Na2SO4"]
        γ_NaNO3(t), [description = "Activity coefficient of NaNO3"]
        γ_NH42SO4(t), [description = "Activity coefficient of (NH4)2SO4"]
        γ_NH4NO3(t), [description = "Activity coefficient of NH4NO3"]
        γ_NH4Cl(t), [description = "Activity coefficient of NH4Cl"]
        γ_NH4HSO4(t), [description = "Activity coefficient of NH4HSO4"]
        γ_NaHSO4(t), [description = "Activity coefficient of NaHSO4"]
        γ_NH43HSO42(t), [description = "Activity coefficient of (NH4)3HSO4_2"]
        γ_H2SO4(t), [description = "Activity coefficient of H2SO4"]
        γ_HHSO4(t), [description = "Activity coefficient of HHSO4"]
        γ_HNO3_aqs(t), [description = "Activity coefficient of HNO3 (salt)"]
        γ_HCl_aqs(t), [description = "Activity coefficient of HCl (salt)"]

        a_CaNO32(t), [description = "Activity of CaNO3_2", unit = u"mol^3/kg^3"]
        a_CaCl2(t), [description = "Activity of CaCl2", unit = u"mol^3/kg^3"]
        a_CaSO4(t), [description = "Activity of CaSO4", unit = u"mol^2/kg^2"]
        a_KHSO4(t), [description = "Activity of KHSO4", unit = u"mol^2/kg^2"]
        a_K2SO4(t), [description = "Activity of K2SO4", unit = u"mol^3/kg^3"]
        a_KNO3(t), [description = "Activity of KNO3", unit = u"mol^2/kg^2"]
        a_KCl(t), [description = "Activity of KCl", unit = u"mol^2/kg^2"]
        a_MgSO4(t), [description = "Activity of MgSO4", unit = u"mol^2/kg^2"]
        a_MgNO32(t), [description = "Activity of MgNO3_2", unit = u"mol^3/kg^3"]
        a_MgCl2(t), [description = "Activity of MgCl2", unit = u"mol^3/kg^3"]
        a_NaCl(t), [description = "Activity of NaCl", unit = u"mol^2/kg^2"]
        a_Na2SO4(t), [description = "Activity of Na2SO4", unit = u"mol^3/kg^3"]
        a_NaNO3(t), [description = "Activity of NaNO3", unit = u"mol^2/kg^2"]
        a_NH42SO4(t), [description = "Activity of (NH4)2SO4", unit = u"mol^2/kg^2"]
        a_NH4NO3(t), [description = "Activity of NH4NO3", unit = u"mol^2/kg^2"]
        a_NH4Cl(t), [description = "Activity of NH4Cl", unit = u"mol^2/kg^2"]
        a_NH4HSO4(t), [description = "Activity of NH4HSO4", unit = u"mol^2/kg^2"]
        a_NaHSO4(t), [description = "Activity of NaHSO4", unit = u"mol^2/kg^2"]
        a_NH43HSO42(t), [description = "Activity of (NH4)3HSO4_2", unit = u"mol^5/kg^5"]
        a_H2SO4(t), [description = "Activity of H2SO4", unit = u"mol^3/kg^3"]
        a_HHSO4(t), [description = "Activity of HHSO4", unit = u"mol^2/kg^2"]
        a_HNO3_aqs(t), [description = "Activity of HNO3 (salt)", unit = u"mol^2/kg^2"]
        a_HCl_aqs(t), [description = "Activity of HCl (salt)", unit = u"mol^2/kg^2"]
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
        HNO3 = Ion(z = 0)
        NH3 = Ion(z = 0)
        HCl = Ion(z = 0)
        HSO4 = Ion(z = abs(-1))
        OH = Ion(z = abs(-1))

        # Salts
        CaNO32 = Salt(cation = Ca, anion = NO3, drh = 0.4906, l_t = 509.4, q = 0.93)
        CaCl2 = Salt(cation = Ca, anion = Cl, drh = 0.2830, l_t = 551.1, q = 2.4)
        CaSO4 = Salt(cation = Ca, anion = SO4, drh = 0.9700, l_t = NaN, q = q0)
        KHSO4 = Salt(cation = K, anion = HSO4, drh = 0.8600, l_t = NaN, q = q0)
        K2SO4 = Salt(cation = K, anion = SO4, drh = 0.9751, l_t = 35.6, q = -0.25)
        KNO3 = Salt(cation = K, anion = NO3, drh = 0.9248, l_t = NaN, q = -2.33)
        KCl = Salt(cation = K, anion = Cl, drh = 0.8426, l_t = 158.9, q = 0.92)
        MgSO4 = Salt(cation = Mg, anion = SO4, drh = 0.8613, l_t = -714.5, q = 0.15)
        MgNO32 = Salt(cation = Mg, anion = NO3, drh = 0.5400, l_t = 230.2, q = 2.32)
        MgCl2 = Salt(cation = Mg, anion = Cl, drh = 0.3284, l_t = 42.23, q = 2.90)
        NaCl = Salt(cation = Na, anion = Cl, drh = 0.7528, l_t = 25.0, q = 2.23)
        Na2SO4 = Salt(cation = Na, anion = SO4, drh = 0.9300, l_t = 80.0, q = -0.19)
        NaNO3 = Salt(cation = Na, anion = NO3, drh = 0.7379, l_t = 304.0, q = -0.39)
        NH42SO4 = Salt(cation = NH4, anion = SO4, drh = 0.7997, l_t = 80.0, q = -0.25)
        NH4NO3 = Salt(cation = NH4, anion = NO3, drh = 0.6183, l_t = 852.0, q = -1.15)
        NH4Cl = Salt(cation = NH4, anion = Cl, drh = 0.7710, l_t = 239.0, q = 0.82)
        NH4HSO4 = Salt(cation = NH4, anion = HSO4, drh = 0.4000, l_t = 384.0, q = q0)
        NaHSO4 = Salt(cation = Na, anion = HSO4, drh = 0.5200, l_t = -45.0, q = q0)
        NH43HSO42 = Salt(cation = NH4, anion = HSO4, drh = 0.6900, l_t = 186.0, q = q0)
        H2SO4 = Salt(cation = H, anion = SO4, drh = 0.000, l_t = NaN, q = -0.1)
        HHSO4 = Salt(cation = H, anion = HSO4, drh = 0.000, l_t = NaN, q = 8.00)
        HNO3_aqs = Salt(cation = H, anion = NO3, drh = NaN, l_t = NaN, q = 2.60) # There is no aqueous to solid conversion for HNO3.
        HCl_aqs = Salt(cation = H, anion = Cl, drh = NaN, l_t = NaN, q = 6.00) # There is no aqueous to solid conversion for HCl.

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
        # Equation 6
        logγₜ₀_CaNO32 ~ -Aᵧ * CaNO32.zz * √I / (√I_one + √I) / √I_one + # NOTE: The last √I_one here is not in the paper but is needed to make the units balance.
                        (CaNO32.zz / (Ca.z + NO3.z)) * (F_Ca / Ca.z + F_NO3 / NO3.z)
        logγₜ₀_CaCl2 ~ -Aᵧ * CaCl2.zz * √I / (√I_one + √I) / √I_one +
                       (CaCl2.zz / (Ca.z + Cl.z)) * (F_Ca / Ca.z + F_Cl / Cl.z)
        logγₜ₀_K2SO4 ~ -Aᵧ * K2SO4.zz * √I / (√I_one + √I) / √I_one +
                       (K2SO4.zz / (2K.z + SO4.z)) * (F_K / K.z + F_SO4 / SO4.z)
        logγₜ₀_KNO3 ~ -Aᵧ * KNO3.zz * √I / (√I_one + √I) / √I_one +
                      (KNO3.zz / (K.z + NO3.z)) * (F_K / K.z + F_NO3 / NO3.z)
        logγₜ₀_KCl ~ -Aᵧ * KCl.zz * √I / (√I_one + √I) / √I_one +
                     (KCl.zz / (K.z + Cl.z)) * (F_K / K.z + F_Cl / Cl.z)
        logγₜ₀_MgSO4 ~ -Aᵧ * MgSO4.zz * √I / (√I_one + √I) / √I_one +
                       (MgSO4.zz / (Mg.z + SO4.z)) * (F_Mg / Mg.z + F_SO4 / SO4.z)
        logγₜ₀_MgNO32 ~ -Aᵧ * MgNO32.zz * √I / (√I_one + √I) / √I_one +
                        (MgNO32.zz / (Mg.z + NO3.z)) * (F_Mg / Mg.z + F_NO3 / NO3.z)
        logγₜ₀_MgCl2 ~ -Aᵧ * MgCl2.zz * √I / (√I_one + √I) / √I_one +
                       (MgCl2.zz / (Mg.z + Cl.z)) * (F_Mg / Mg.z + F_Cl / Cl.z)
        logγₜ₀_NaCl ~ -Aᵧ * NaCl.zz * √I / (√I_one + √I) / √I_one +
                      (NaCl.zz / (Na.z + Cl.z)) * (F_Na / Na.z + F_Cl / Cl.z)
        logγₜ₀_Na2SO4 ~ -Aᵧ * Na2SO4.zz * √I / (√I_one + √I) / √I_one +
                        (Na2SO4.zz / (2Na.z + SO4.z)) * (F_Na / Na.z + F_SO4 / SO4.z)
        logγₜ₀_NaNO3 ~ -Aᵧ * NaNO3.zz * √I / (√I_one + √I) / √I_one +
                       (NaNO3.zz / (Na.z + NO3.z)) * (F_Na / Na.z + F_NO3 / NO3.z)
        logγₜ₀_NH42SO4 ~ -Aᵧ * NH42SO4.zz * √I / (√I_one + √I) / √I_one +
                         (NH42SO4.zz / (2NH4.z + SO4.z)) * (F_NH4 / NH4.z + F_SO4 / SO4.z)
        logγₜ₀_NH4NO3 ~ -Aᵧ * NH4NO3.zz * √I / (√I_one + √I) / √I_one +
                        (NH4NO3.zz / (NH4.z + NO3.z)) * (F_NH4 / NH4.z + F_NO3 / NO3.z)
        logγₜ₀_NH4Cl ~ -Aᵧ * NH4Cl.zz * √I / (√I_one + √I) / √I_one +
                       (NH4Cl.zz / (NH4.z + Cl.z)) * (F_NH4 / NH4.z + F_Cl / Cl.z)
        logγₜ₀_H2SO4 ~ -Aᵧ * H2SO4.zz * √I / (√I_one + √I) / √I_one +
                       (H2SO4.zz / (2H.z + SO4.z)) * (F_H / H.z + F_SO4 / SO4.z)
        logγₜ₀_HHSO4 ~ -Aᵧ * HHSO4.zz * √I / (√I_one + √I) / √I_one +
                       (HHSO4.zz / (H.z + HSO4.z)) * (F_H / H.z + F_HSO4 / HSO4.z)
        logγₜ₀_HNO3_aqs ~ -Aᵧ * HNO3_aqs.zz * √I / (√I_one + √I) / √I_one +
                          (HNO3_aqs.zz / (H.z + NO3.z)) * (F_H / H.z + F_NO3 / NO3.z)
        logγₜ₀_HCl_aqs ~ -Aᵧ * HCl_aqs.zz * √I / (√I_one + √I) / √I_one +
                         (HCl_aqs.zz / (H.z + Cl.z)) * (F_H / H.z + F_Cl / Cl.z)

        Aᵧ_term ~ Aᵧ * √I / (√I_one + √I) / √I_one # NOTE: The last √I_one here is not in the paper but is needed to make the units balance.
        # Equation 7
        F_NH4 ~ sum([s.Y * s.logγ⁰ for s in [NH4NO3, NH4Cl, NH4HSO4, NH42SO4, NH43HSO42]]) +
                Aᵧ_term *
                sum([s.zz * s.Y for s in [NH4NO3, NH4Cl, NH4HSO4, NH42SO4, NH43HSO42]])
        F_Na ~ sum([s.Y * s.logγ⁰ for s in [NaCl, Na2SO4, NaNO3, NaHSO4]]) +
               Aᵧ_term * sum([s.zz * s.Y for s in [NaCl, Na2SO4, NaNO3, NaHSO4]])
        F_H ~ sum([s.Y * s.logγ⁰ for s in [H2SO4, HHSO4, HCl_aqs, HNO3_aqs]]) +
              Aᵧ_term * sum([s.zz * s.Y for s in [H2SO4, HHSO4, HCl_aqs, HNO3_aqs]])
        F_Ca ~ sum([s.Y * s.logγ⁰ for s in [CaNO32, CaCl2, CaSO4]]) +
               Aᵧ_term * sum([s.zz * s.Y for s in [CaNO32, CaCl2, CaSO4]])
        F_K ~ sum([s.Y * s.logγ⁰ for s in [KHSO4, K2SO4, KNO3, KCl]]) +
              Aᵧ_term * sum([s.zz * s.Y for s in [KHSO4, K2SO4, KNO3, KCl]])
        F_Mg ~ sum([s.Y * s.logγ⁰ for s in [MgSO4, MgNO32, MgCl2]]) +
               Aᵧ_term * sum([s.zz * s.Y for s in [MgSO4, MgNO32, MgCl2]])

        # Equation 8
        F_Cl ~ sum([s.X * s.logγ⁰ for s in [NaCl, KCl, MgCl2, CaCl2, NH4Cl, HCl_aqs]]) +
               Aᵧ_term *
               sum([s.zz * s.X for s in [NaCl, KCl, MgCl2, CaCl2, NH4Cl, HCl_aqs]])
        F_NO3 ~ sum([s.X * s.logγ⁰ for s in [NaNO3, KNO3, MgNO32, CaNO32, NH4NO3, HNO3_aqs]]) +
                Aᵧ_term *
                sum([s.zz * s.X for s in [NaNO3, KNO3, MgNO32, CaNO32, NH4NO3, HNO3_aqs]])
        F_SO4 ~ sum([s.X * s.logγ⁰ for s in [Na2SO4, K2SO4, MgSO4, CaSO4, NH42SO4, H2SO4]]) +
                Aᵧ_term *
                sum([s.zz * s.X for s in [Na2SO4, K2SO4, MgSO4, CaSO4, NH42SO4, H2SO4]])
        F_HSO4 ~ sum([s.X * s.logγ⁰ for s in [KHSO4, HHSO4, NaHSO4, NH4HSO4, NH43HSO42]]) +
                 Aᵧ_term *
                 sum([s.zz * s.X for s in [KHSO4, HHSO4, NaHSO4, NH4HSO4, NH43HSO42]])

        # Equation 14
        logγ_CaNO32 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_CaNO32 -
                      (0.125 - c_1 * (T - T₀₂)) * A
        logγ_CaCl2 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_CaCl2 -
                     (0.125 - c_1 * (T - T₀₂)) * A
        logγ_K2SO4 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_K2SO4 -
                     (0.125 - c_1 * (T - T₀₂)) * A
        logγ_KNO3 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_KNO3 -
                    (0.125 - c_1 * (T - T₀₂)) * A
        logγ_KCl ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_KCl -
                   (0.125 - c_1 * (T - T₀₂)) * A
        logγ_MgSO4 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_MgSO4 -
                     (0.125 - c_1 * (T - T₀₂)) * A
        logγ_MgNO32 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_MgNO32 -
                      (0.125 - c_1 * (T - T₀₂)) * A
        logγ_MgCl2 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_MgCl2 -
                     (0.125 - c_1 * (T - T₀₂)) * A
        logγ_NaCl ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_NaCl -
                    (0.125 - c_1 * (T - T₀₂)) * A
        logγ_Na2SO4 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_Na2SO4 -
                      (0.125 - c_1 * (T - T₀₂)) * A
        logγ_NaNO3 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_NaNO3 -
                     (0.125 - c_1 * (T - T₀₂)) * A
        logγ_NH42SO4 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_NH42SO4 -
                       (0.125 - c_1 * (T - T₀₂)) * A
        logγ_NH4NO3 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_NH4NO3 -
                      (0.125 - c_1 * (T - T₀₂)) * A
        logγ_NH4Cl ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_NH4Cl -
                     (0.125 - c_1 * (T - T₀₂)) * A
        logγ_H2SO4 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_H2SO4 -
                     (0.125 - c_1 * (T - T₀₂)) * A
        logγ_HHSO4 ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_HHSO4 -
                     (0.125 - c_1 * (T - T₀₂)) * A
        logγ_HNO3_aqs ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_HNO3_aqs -
                        (0.125 - c_1 * (T - T₀₂)) * A
        logγ_HCl_aqs ~ (1.125 - c_1 * (T - T₀₂)) * logγₜ₀_HCl_aqs -
                       (0.125 - c_1 * (T - T₀₂)) * A

        # Equation in text below Equation 14
        A ~ -((0.41√I / (√I_one + √I)) + 0.039(I / I_one)^0.92)

        # Activity coefficients
        γ_CaNO32 ~ exp(logγ_CaNO32)
        γ_CaCl2 ~ exp(logγ_CaCl2)
        γ_K2SO4 ~ exp(logγ_K2SO4)
        γ_KNO3 ~ exp(logγ_KNO3)
        γ_KCl ~ exp(logγ_KCl)
        γ_MgSO4 ~ exp(logγ_MgSO4)
        γ_MgNO32 ~ exp(logγ_MgNO32)
        γ_MgCl2 ~ exp(logγ_MgCl2)
        γ_NaCl ~ exp(logγ_NaCl)
        γ_Na2SO4 ~ exp(logγ_Na2SO4)
        γ_NaNO3 ~ exp(logγ_NaNO3)
        γ_NH42SO4 ~ exp(logγ_NH42SO4)
        γ_NH4NO3 ~ exp(logγ_NH4NO3)
        γ_NH4Cl ~ exp(logγ_NH4Cl)
        γ_H2SO4 ~ exp(logγ_H2SO4)
        γ_HHSO4 ~ exp(logγ_HHSO4)
        γ_HNO3_aqs ~ exp(logγ_HNO3_aqs)
        γ_HCl_aqs ~ exp(logγ_HCl_aqs)

        # Table 4 footnotes
        γ_CaSO4 ~ 0.0 # From Table 4 footnote a, CaSO4 has an activity coefficient of zero.
        γ_KHSO4 ~ √(γ_HHSO4 * γ_KCl / γ_HCl_aqs) # From Table 4 footnote b.
        γ_NH4HSO4 ~ √(γ_HHSO4 * γ_NH4Cl / γ_HCl_aqs) # From Table 4 footnote c
        γ_NaHSO4 ~ √(γ_HHSO4 * γ_NaCl / γ_HCl_aqs) # From Table 4 footnote d
        γ_NH43HSO42 ~ (γ_NH42SO4^3 * γ_NH4HSO4)^(1 / 5) # From Table 4 footnote e

        # Ionic strength (equation in text below Equation 8)
        I ~ 0.5 * sum([ion.m * ion.z^2
                 for ion in [
            NH4, Na, H, Ca, K, Mg, Cl, NO3, SO4, HNO3, NH3, HCl, HSO4, OH]])
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
        HNO3_aqs.I ~ I
        HCl_aqs.I ~ I

        # Activities (Section 2.2)
        a_CaNO32 ~ Ca.m * NO3.m^2 * γ_CaNO32^(1 + 2)
        a_CaCl2 ~ Ca.m * Cl.m^2 * γ_CaCl2^(1 + 2)
        a_CaSO4 ~ Ca.m * SO4.m * γ_CaSO4^(1 + 1)
        a_KHSO4 ~ K.m * HSO4.m * γ_KHSO4^(1 + 1)
        a_K2SO4 ~ K.m^2 * SO4.m * γ_K2SO4^(2 + 1)
        a_KNO3 ~ K.m * NO3.m * γ_KNO3^(1 + 1)
        a_KCl ~ K.m * Cl.m * γ_KCl^(1 + 1)
        a_MgSO4 ~ Mg.m * SO4.m * γ_MgSO4^(1 + 1)
        a_MgNO32 ~ Mg.m * NO3.m^2 * γ_MgNO32^(1 + 2)
        a_MgCl2 ~ Mg.m * Cl.m^2 * γ_MgCl2^(1 + 2)
        a_NaCl ~ Na.m * Cl.m * γ_NaCl^(1 + 1)
        a_Na2SO4 ~ Na.m^2 * SO4.m * γ_Na2SO4^(2 + 1)
        a_NaNO3 ~ Na.m * NO3.m * γ_NaNO3^(1 + 1)
        a_NH42SO4 ~ NH4.m * SO4.m * γ_NH42SO4^(1 + 1)
        a_NH4NO3 ~ NH4.m * NO3.m * γ_NH4NO3^(1 + 1)
        a_NH4Cl ~ NH4.m * Cl.m^1 * γ_NH4Cl^(1 + 1)
        a_NH4HSO4 ~ NH4.m * HSO4.m * γ_NH4HSO4^(1 + 1)
        a_NaHSO4 ~ Na.m * HSO4.m * γ_NaHSO4^(1 + 1)
        a_NH43HSO42 ~ NH4.m^3 * HSO4.m^2 * γ_NH43HSO42^(3 + 2)
        a_H2SO4 ~ H.m^2 * SO4.m * γ_H2SO4^(2 + 1)
        a_HHSO4 ~ H.m * HSO4.m * γ_HHSO4^(1 + 1)
        a_HNO3_aqs ~ H.m * NO3.m * γ_HNO3_aqs^(1 + 1)
        a_HCl_aqs ~ H.m * Cl.m * γ_HCl_aqs^(1 + 1)

        # Water content (Section 2.3)
        W ~ CaNO32.m * W / maw_CaNO32.m_aw + CaCl2.m * W / maw_CaCl2.m_aw + KHSO4.m * W / maw_KHSO4.m_aw +
            K2SO4.m * W / maw_K2SO4.m_aw + KNO3.m * W / maw_KNO3.m_aw + KCl.m * W / maw_KCl.m_aw +
            MgSO4.m * W / maw_MgSO4.m_aw + MgNO32.m * W / maw_MgNO32.m_aw + MgCl2.m * W / maw_MgCl2.m_aw +
            NaNO3.m * W / maw_NaNO3.m_aw + NaHSO4.m * W / maw_NaHSO4.m_aw + NaCl.m * W / maw_NaCl.m_aw +
            Na2SO4.m * W / maw_Na2SO4.m_aw + NH42SO4.m * W / maw_NH42SO4.m_aw + NH4NO3.m * W / maw_NH4NO3.m_aw +
            NH4Cl.m * W / maw_NH4Cl.m_aw + NH4HSO4.m * W / maw_NH4HSO4.m_aw + NH43HSO42.m * W / maw_NH43HSO42.m_aw
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

        # Mass balance (need to double check)
        Ca.m ~ CaNO32.m + CaCl2.m + CaSO4.m
        K.m ~ KHSO4.m + K2SO4.m + KNO3.m + KCl.m
        Mg.m ~ MgSO4.m + MgNO32.m + MgCl2.m
        Na.m ~ NaCl.m + Na2SO4.m + NaNO3.m + NaHSO4.m
        NH4.m ~ NH4NO3.m + NH4Cl.m + NH4HSO4.m + 2NH42SO4.m + 3NH43HSO42.m
        H.m ~ 2H2SO4.m + HHSO4.m + HCl_aqs.m + HNO3_aqs.m
        Cl.m ~ NaCl.m + KCl.m + MgCl2.m + CaCl2.m + NH4Cl.m + HCl_aqs.m
        NO3.m ~ NaNO3.m + KNO3.m + MgNO32.m + CaNO32.m + NH4NO3.m + HNO3_aqs.m
        SO4.m ~ Na2SO4.m + K2SO4.m + MgSO4.m + CaSO4.m + NH42SO4.m + H2SO4.m
        HSO4.m ~ KHSO4.m + HHSO4.m + NaHSO4.m + NH4HSO4.m + NH43HSO42.m

        # Charge balance
    end
end

@named aq = Aqueous()
equations(aq)

@mtkmodel AqueousTest begin
    @components begin
        aq = Aqueous()
    end
    @constants begin
        no_change = 0.0, [unit = u"mol/kg/s"]
    end
    @equations begin
        D(aq.Ca.m) ~ no_change
        D(aq.Cl.m) ~ no_change
        D(aq.SO4.m) ~ no_change
        D(aq.Na.m) ~ no_change
        D(aq.NH4.m) ~ no_change
        D(aq.H.m) ~ no_change
        D(aq.K.m) ~ no_change
        D(aq.Mg.m) ~ no_change
        D(aq.NO3.m) ~ no_change
        D(aq.HSO4.m) ~ no_change
        D(aq.OH.m) ~ no_change
        D(aq.HNO3.m) ~ no_change
        D(aq.NH3.m) ~ no_change
        D(aq.HCl.m) ~ no_change
    end
end

@named aqt = AqueousTest()

sys = mtkcompile(aqt)
unknowns(sys)

prob = ODEProblem(sys, [], (0.0, 1.0))

using SymbolicIndexingInterface: setp, getsym, parameter_values
using SciMLBase: remake
prob = remake(prob, u0 = [sys.aq.Ca.m => 1.0])
f = getsym(prob, [sys.aq.γ_NaCl, sys.aq.γ_CaCl2, sys.aq.γ_NaNO3, sys.aq.γ_CaNO32])
f(prob)

prob = remake(prob, u0 = [sys.aq.Cl.m => 1.0])
f = getsym(prob, [sys.aq.γ_NaCl, sys.aq.γ_CaCl2, sys.aq.γ_NaNO3, sys.aq.γ_CaNO32])
f(prob)
