@constants tinyrate = 1.e-99 [unit = u"mol/m_air^3/s", description = "Tiny reaction rate constant negative concentrations"]

"""
Define a reaction based on information in Table 2 of Fountoukis and Nenes (2007).
"""
struct Rxn
    reactant
    product
    sys

    function Rxn(reactant, product, K⁰::Number, K⁰units, hgroup::Number, cgroup::Number; name="rxn")
        γr = γ(reactant)
        γp = γ(product)
        ar = activity(reactant)
        ap = activity(product)
        # These are the variables from Fountoukis and Nenes (2007) Table 2
        @constants K⁰ = K⁰ [unit = K⁰units, description = "Equilibrium constant at 298.15 K"] 
        @constants H_group = hgroup [description = "ΔH⁰ / (R * T₀) (unitless)"]
        @constants C_group = cgroup [description = "ΔC⁰ₚ / R (unitless)"]
        @variables K_eq(t) [unit = K⁰units, description = "Equilibrium constant"]
        # These are the transformed variables to turn this into a kinetic reaction rather than an equilibrium reaction.
        # To do so we need to choose a base reaction rate constant to make reactions proceed at a reasonable rate.
        baserate = 1
        @constants fakerate = baserate [unit = u"mol/m_air^3/s", description = "Fake reaction rate constant to convert equilibrium reaction into kinetic reaction"]
        @variables γ_r(t) [unit=ModelingToolkit.get_unit(γr), description="Activity coefficient of reactant"]
        @variables γ_p(t) [unit=ModelingToolkit.get_unit(γp), description="Activity coefficient of product"]
        @constants k_fwd=baserate [unit = ModelingToolkit.get_unit(γr/ar*fakerate), description = "Forward reaction rate constant"]
        @variables k_rev(t)=1 [unit = ModelingToolkit.get_unit(γp/ap*fakerate), description = "Reverse reaction rate constant"]
        @constants unit_krev=1 [unit = ModelingToolkit.get_unit(γp/ap*fakerate), description = "Unit reverse reaction rate constant"]
        if (typeof(product) <: SaltLike) && (reactant isa Solid) # Deliquescence
            krev = f_deliquescence[product] * γ_p / γ_r / K_eq * fakerate
        else # Not an aqueous precipitation reaction so no deliquescence
            krev = γ_p / γ_r / K_eq * fakerate
        end
        @constants unitconv = 1 [unit = ModelingToolkit.get_unit(k_rev/krev), description = "Unit conversion for k_rev"]
        sys = NonlinearSystem([
            K_eq ~ K⁰ * exp(-H_group * (T₀ / T - 1) - C_group * (1 + log(T₀ / T) - T₀ / T))
            γ_r ~ γr
            γ_p ~ γp
            k_rev ~ min(unit_krev*1e20, krev * unitconv)
        ], [K_eq, γ_r, γ_p, k_rev], [K⁰, H_group, C_group, k_fwd, fakerate, unitconv, unit_krev, T₀, T₀₂, c_1, I_one]; name=name)
        new(reactant, product, sys)
    end
end

# Equation 5: Equilibrium constant
@constants T₀ = 293.15 [unit = u"K", description = "Standard temperature"]

"""
Assemble an equation for a reaction based on Table 2 of Fountoukis and Nenes (2007), where
the left-hand side is the equilibrium constant and the right-hand side is the activity.
"""
function rxn_eqs(r::Rxn)
    rterms = terms(r.reactant)
    pterms = terms(r.product)
    fwd = Reaction(r.sys.k_fwd,  rterms[1], pterms[1], rterms[2], pterms[2])
    rev = Reaction(r.sys.k_rev, pterms[1], rterms[1], pterms[2], rterms[2])
    return [fwd, rev]
end

# NOTE: Assuming that H_group and C_group are zero when they are left out of Table 2. 
@named rxn1 = Rxn(CaNO32s, CaNO32_aqs, 6.067e5, u"mol^3/kg_water^3", -11.299, 0.0)
@named rxn2 = Rxn(CaCl2s, CaCl2_aqs, 7.974e11, u"mol^3/kg_water^3", -14.087, 0.0)
@named rxn3 = Rxn(CaSO4s, CaSO4_aqs, 4.319e-5, u"mol^2/kg_water^2", 0.0, 0.0)
@named rxn4 = Rxn(K2SO4s, K2SO4_aqs, 1.569e-2, u"mol^3/kg_water^3", -9.589, 45.807)
@named rxn5 = Rxn(KHSO4s, KHSO4_aqs, 24.016, u"mol^2/kg_water^2", -8.423, 17.964)
@named rxn6 = Rxn(KNO3s, KNO3_aqs, 0.872, u"mol^2/kg_water^2", 14.075, 19.388)
@named rxn7 = Rxn(KCls, KCl_aqs, 8.680, u"mol^2/kg_water^2", -6.167, 19.953)
@named rxn8 = Rxn(MgSO4s, MgSO4_aqs, 1.079e5, u"mol^2/kg_water^2", 36.798, 0.0)
@named rxn9 = Rxn(MgNO32s, MgNO32_aqs, 2.507e15, u"mol^3/kg_water^3", -8.754, 0.0)
@named rxn10 = Rxn(MgCl2s, MgCl2_aqs, 9.557e21, u"mol^3/kg_water^3", -1.347, 0.0)
@named rxn11 = Rxn(HSO4_ion, [H_ion, SO4_ion], 1.015e-2, u"mol/kg_water", 8.85, 25.14)
@named rxn12 = Rxn(NH3g, NH3_ion, 5.764e1, u"mol/kg_water/atm", 13.79, -5.39)
@named rxn13 = Rxn([NH3_ion, H2Oaq], [NH4_ion, OH_ion], 1.805e-5, u"mol/kg_water", -1.50, 26.92)
@named rxn14 = Rxn(HNO3g, HNO3_aqs, 2.511e6, u"mol^2/kg_water^2/atm", 29.17, 16.83)
@named rxn15 = Rxn(HNO3g, HNO3_ion, 2.1e5, u"mol/kg_water/atm", 29.17, 16.83)
@named rxn16 = Rxn(HClg, [H_ion, Cl_ion], 1.971e6, u"mol^2/kg_water^2/atm", 30.20, 19.91)
@named rxn17 = Rxn(HClg, HCl_ion, 2.5e3, u"mol/kg_water/atm", 30.20, 19.91)
@named rxn18 = Rxn(H2Oaq, [H_ion, OH_ion], 1.010e-14, u"mol^2/kg_water^2", -22.52, 26.92)
@named rxn19 = Rxn(Na2SO4s, Na2SO4_aqs, 4.799e-1, u"mol^3/kg_water^3", 0.98, 39.75)
@named rxn20 = Rxn(NH42SO4s, NH42SO4_aqs, 1.87e0, u"mol^3/kg_water^3", -2.65, 38.57)
@named rxn21 = Rxn(NH4Cls, [NH3g, HClg], 1.086e-16, u"atm^2", -71.00, 2.40)
@named rxn22 = Rxn(NaNO3s, NaNO3_aqs, 1.197e1, u"mol^2/kg_water^2", -8.22, 16.01)
@named rxn23 = Rxn(NaCls, NaCl_aqs, 3.766e1, u"mol^2/kg_water^2", -1.56, 16.90)
@named rxn24 = Rxn(NaHSO4s, NaHSO4_aqs, 2.413e4, u"mol^2/kg_water^2", 0.79, 14.75)
@named rxn25 = Rxn(NH4NO3s, [NH3g, HNO3g], 4.199e-17, u"atm^2", -74.375, 6.025)
@named rxn26 = Rxn(NH4HSO4s, NH4HSO4_aqs, 1.383e0, u"mol^2/kg_water^2", -2.87, 15.83)
@named rxn27 = Rxn(NH43HSO42s, NH43HSO42_aqs, 2.972e1, u"mol^5/kg_water^5", -5.19, 54.40)

all_rxns = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6, rxn7, rxn8, rxn9, rxn10, rxn11, rxn12, rxn13, rxn14,
    rxn15, rxn16, rxn17, rxn18, rxn19, rxn20, rxn21, rxn22, rxn23, rxn24, rxn25, rxn26, rxn27]

all_rxn_eqs = [rxn_eqs(x) for x in all_rxns]

@named IW = NonlinearSystem([     
        # Calculate the ionic strength of the multicomponent solution as described by 
        # Fountoukis and Nenes (2007), between equations 8 and 9: ``I = \\frac{1}{2} \\sum_i m_i z_i^2``
        # Force I to always be positive to avoid attempts to take the square root of a negative number.
        I ~ max(1.0e-20*I_one, 1 / 2 * sum([ion.m * ion.z^2 for ion in all_Ions] / W))

        # Water content (Fountoukis and Nenes Eq. 16).
        W ~ max(1.0e-10*W_one, W_eq16)
    ], [I; W], [])