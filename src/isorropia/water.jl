struct H2O <: Species 
    RH
end

Base.nameof(w::H2O) = "H2O"

vars(w::H2O) = []
terms(w::H2O) = [], []

# From Equation 15 in Fountoukis and Nenes (2007), the activity of water is 
# equal to the relative humidity
activity(w::H2O) = w.RH

@constants unit_molality=1.0 [unit = u"mol/kg_water", description = "Unit molality"]

"""
Create an equation system for the water content of the aerosol based on Fountoukis and Nenes (2007) Eq. 16.
"""
function Water(t, active_salts)
    # Polynomial fit for molality at given water activity a_w from Table 7 of Fountoukis and Nenes (2007).
    m_aw_coeffs = Dict(
        :CaNO32_aqs => [36.356, -165.66, 447.46, -673.55, 510.91, -155.56, 0] * unit_molality,
        :CaCl2_aqs => [20.847, -97.599, 273.220, -422.120, 331.160, -105.450, 0] * unit_molality,
        :KHSO4_aqs => [1.061, -0.101, 1.579 * 10^-2, -1.950 * 10^-3, 9.515 * 10^-5, -1.547 * 10^-6, 0] * unit_molality,
        :K2SO4_aqs => [1061.51, -4748.97, 8096.16, -6166.16, 1757.47, 0, 0] * unit_molality,
        :KNO3_aqs => [1.2141 * 10^4, -5.1173 * 10^4, 8.1252 * 10^4, -5.7527 * 10^4, 1.5305 * 10^4, 0, 0] * unit_molality,
        :KCl_aqs => [179.721, -721.266, 1161.03, -841.479, 221 / 943, 0, 0] * unit_molality,
        :MgSO4_aqs => [-0.778, 177.74, -719.79, 1174.6, -863.44, 232.31, 0] * unit_molality,
        :MgNO32_aqs => [12.166, -16.154, 0, 10.886, 0, -6.815, 0] * unit_molality,
        :MgCl2_aqs => [11.505, -26.518, 34.937, -19.829, 0, 0, 0] * unit_molality,
        :NaNO3_aqs => [0.9988, -2.6947*10^-2, 1.9610*10^-4, 2.8154*10^-5, 6.1359*10^-7, 0, 0] * unit_molality,
        :NaHSO4_aqs => [1.0614, -0.1014, 1.5796*10^-2, -1.9501*10^-3, 9.5147*10^-5, -1.5473*10^-6, 0] * unit_molality,
        :NaCl_aqs => [1.0084, -4.9390*10^-2, 8.888*10^-3, -2.1570*10^-3, 1.6170*10^-4, 1.99*10^-6, -1.142*10^-7] * unit_molality,
        :Na2SO4_aqs => [1.0052, -6.4840*10^-2, 3.519*10^-2, -1.3190*10^-2, 1.9250*10^-3, -1.224*10^-4, 2.87*10^-6] * unit_molality,
        :NH42SO4_aqs => [0.9968, -2.9690*10^-2, 1.735*10^-5, -3.2530*10^-4, 3.5710*10^-5, -9.7870*10^-7, 0] * unit_molality,
        :NH4Cl_aqs => [0.9968, -2.6110*10^-2, -1.5990*10^-3, 1.3550*10^-4, -2.317*10^-6, -1.113*10^-8, 0] * unit_molality,
        :NH4NO3_aqs => [1.0053, -2.4991*10^-2, 4.4688*10^-4, 1.6453*10^-5, -3.8940*10^-7, -4.7668*10^-8, 1.3753*10^-9] * unit_molality,
        :NH4HSO4_aqs => [1.0261, -4.9766*10^-2, 3.2757*10^-3, -2.4477*10^-4, 1.0766*10^-5, -1.8329*10^-7, 0] * unit_molality,
        :NH43HSO42_aqs => [1.0088, -5.3730*10^-2, 1.4201*10^-3, -9.2484*10^-4, 2.2796*10^-4, -1.5445*10^-5, 0] * unit_molality
    )
    @variables W(t) = 1.0e-8 [unit = u"kg_water/m_air^3", description = "Aerosol water content"]
    W = ParentScope(W)
    w = 0
    for s ∈ active_salts
        n = Symbol(nameof(s))
        if n ∈ keys(m_aw_coeffs)
            m_aw_coeff = m_aw_coeffs[n]
            # Water content (Fountoukis and Nenes Eq. 16).
            w += (s.cation.m + s.anion.m) / sum(m_aw_coeff .* RH.^collect(0:6)) 
        end
    end
    @named water = ODESystem([W ~ w], t, [W], [])
end