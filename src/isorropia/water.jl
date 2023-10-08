struct H2O <: Species 
    m
end

# From Equation 15 in Fountoukis and Nenes (2007), the activity of water is 
# equal to the relative humidity
@parameters H2O_aq(t) = 1e-9 [unit = u"mol/m_air^3", isconstantspecies=true,
        description = "Dummy water concentration to make units balance (the real water concentration is variable `W`)"]
H2O_aq = ParentScope(H2O_aq)
H2Oaq = H2O(H2O_aq)

γ(w::H2O) = RH / w.m
terms(w::H2O) = [w.m], [1]
min_conc(w::H2O) = unit_conc

@constants unit_molality=1.0 [unit = u"mol/kg_water", description = "Unit molality"]

# Polynomial fit for molality at given water activity a_w from Table 7 of Fountoukis and Nenes (2007).
m_aw_coeffs = [
    CaNO32_aqs => [36.356, -165.66, 447.46, -673.55, 510.91, -155.56, 0] * unit_molality,
    CaCl2_aqs => [20.847, -97.599, 273.220, -422.120, 331.160, -105.450, 0] * unit_molality,
    KHSO4_aqs => [1.061, -0.101, 1.579 * 10^-2, -1.950 * 10^-3, 9.515 * 10^-5, -1.547 * 10^-6, 0] * unit_molality,
    K2SO4_aqs => [1061.51, -4748.97, 8096.16, -6166.16, 1757.47, 0, 0] * unit_molality,
    KNO3_aqs => [1.2141 * 10^4, -5.1173 * 10^4, 8.1252 * 10^4, -5.7527 * 10^4, 1.5305 * 10^4, 0, 0] * unit_molality,
    KCl_aqs => [179.721, -721.266, 1161.03, -841.479, 221 / 943, 0, 0] * unit_molality,
    MgSO4_aqs => [-0.778, 177.74, -719.79, 1174.6, -863.44, 232.31, 0] * unit_molality,
    MgNO32_aqs => [12.166, -16.154, 0, 10.886, 0, -6.815, 0] * unit_molality,
    MgCl2_aqs => [11.505, -26.518, 34.937, -19.829, 0, 0, 0] * unit_molality,
    # TODO(CT): We are missing some salts here.
]

# Equation 16.
W_eq16 = sum([(salt.cation.m / salt.ν_cation) / sum(m_aw_coeff .* RH.^collect(0:6)) 
    for (salt, m_aw_coeff) ∈ m_aw_coeffs])