@mtkmodel EqConst begin
    @description """An equilibrium constant based on Equation 5 in Fountoukis and Nenes (2007),
parameterized by values from Table 2 in the same paper."""
    @structural_parameters begin
        K⁰units = u"mol^3/kg^3"
    end
    @constants begin
        K⁰, [description = "Equilibrium constant at 298.15 K"]
        H⁺, [description = "ΔH⁰ / (R * T₀) (unitless)"]
        C⁺, [description = "ΔC⁰ₚ / R (unitless)"]
        T₀ = 293.15, [unit = u"K", description = "Standard temperature"]
    end
    @parameters begin
        T, [unit = u"K", description = "Temperature"]
    end
    @variables begin
        K_eq(t), [description = "Equilibrium constant"]
    end
    @equations begin
         K_eq ~ K⁰ * exp(-H⁺ * (T₀ / T - 1) - C⁺ * (1 + log(T₀ / T) - T₀ / T))
    end
end

@mtkmodel EquilibriumConstants begin
    @parameters begin
        T, [unit = u"K", description = "Temperature"]
    end
    @components begin
     # Equilibrium constants from Table 2 of Fountoukis and Nenes (2007).
     # NOTE: Assuming that H⁺ and C⁺ are zero when they are left out of Table 2.
     r1 = EqConst(K⁰ = 6.067e5, H⁺ = -11.299, C⁺ = 0.0)
     r2 = EqConst(K⁰ = 7.974e11, H⁺ = -14.087, C⁺ = 0.0)
     r3 = EqConst(K⁰ = 4.319e-5, H⁺ = 0.0, C⁺ = 0.0)
     r4 = EqConst(K⁰ = 1.569e-2, H⁺ = -9.589, C⁺ = 45.807)
     r5 = EqConst(K⁰ = 24.016, H⁺ = -8.423, C⁺ = 17.964)
     r6 = EqConst(K⁰ = 0.872, H⁺ = 14.075, C⁺ = 19.388)
     r7 = EqConst(K⁰ = 8.680, H⁺ = -6.167, C⁺ = 19.953)
     r8 = EqConst(K⁰ = 1.079e5, H⁺ = 36.798, C⁺ = 0.0)
     r9 = EqConst(K⁰ = 2.507e15, H⁺ = -8.754, C⁺ = 0.0)
     r10 = EqConst(K⁰ = 9.557e21, H⁺ = -1.347, C⁺ = 0.0)
     r11 = EqConst(K⁰ = 1.015e-2, H⁺ = 8.85, C⁺ = 25.14)
     r12 = EqConst(K⁰ = 5.764e1, H⁺ = 13.79, C⁺ = -5.39)
     r13 = EqConst(K⁰ = 1.805e-5, H⁺ = -1.50, C⁺ = 26.92)
     r14 = EqConst(K⁰ = 2.511e6, H⁺ = 29.17, C⁺ = 16.83)
     r15 = EqConst(K⁰ = 2.1e5, H⁺ = 29.17, C⁺ = 16.83)
     r16 = EqConst(K⁰ = 1.971e6, H⁺ = 30.20, C⁺ = 19.91)
     r17 = EqConst(K⁰ = 2.5e3, H⁺ = 30.20, C⁺ = 19.91)
     r18 = EqConst(K⁰ = 1.010e-14, H⁺ = -22.52, C⁺ = 26.92)
     r19 = EqConst(K⁰ = 4.799e-1, H⁺ = 0.98, C⁺ = 39.75)
     r20 = EqConst(K⁰ = 1.87e0, H⁺ = -2.65, C⁺ = 38.57)
     r21 = EqConst(K⁰ = 1.086e-16, H⁺ = -71.00, C⁺ = 2.40)
     r22 = EqConst(K⁰ = 1.197e1, H⁺ = -8.22, C⁺ = 16.01)
     r23 = EqConst(K⁰ = 3.766e1, H⁺ = -1.56, C⁺ = 16.90)
     r24 = EqConst(K⁰ = 2.413e4, H⁺ = 0.79, C⁺ = 14.75)
     r25 = EqConst(K⁰ = 4.199e-17, H⁺ = -74.375, C⁺ = 6.025)
     r26 = EqConst(K⁰ = 1.383e0, H⁺ = -2.87, C⁺ = 15.83)
     r27 = EqConst(K⁰ = 2.972e1, H⁺ = -5.19, C⁺ = 54.40)
    end
    @constants begin
        k1_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k2_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k3_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k4_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k5_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k6_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k7_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k8_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k9_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k10_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k11_unit = 1, [unit = u"mol/kg", description = "Unit conversion factor"]
        k12_unit = 1, [unit = u"mol/kg/Constants.atm", description = "Unit conversion"]
        k13_unit = 1, [unit = u"mol/kg", description = "Unit conversion factor"]
        k14_unit = 1, [unit = u"mol^2/kg^2/Constants.atm", description = "Unit conversion"]
        k1b_unit = 1, [unit = u"mol/kg", description = "Unit conversion"]
        k15_unit = 1, [unit = u"mol^2/kg^2/Constants.atm", description = "Unit conversion"]
        k16_unit = 1, [unit = u"mol^2/kg^2/Constants.atm", description = "Unit conversion"]
        k17_unit = 1, [unit = u"mol/kg/Constants.atm", description = "Unit conversion"]
        k2b_unit = 1, [unit = u"mol/kg", description = "Unit conversion factor"]
        k18_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k19_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k20_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k21_unit = 1, [unit = u"Constants.atm^2", description = "Unit conversion factor"]
        k22_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k23_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k24_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k25_unit = 1, [unit = u"Constants.atm^2", description = "Unit conversion factor"]
        k26_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k27_unit = 1, [unit = u"mol^5/kg^5", description = "Unit conversion factor"]
    end
    @equations begin
        r1.T ~ T
        r2.T ~ T
        r3.T ~ T
        r4.T ~ T
        r5.T ~ T
        r6.T ~ T
        r7.T ~ T
        r8.T ~ T
        r9.T ~ T
        r10.T ~ T
        r11.T ~ T
        r12.T ~ T
        r13.T ~ T
        r14.T ~ T
        r15.T ~ T
        r16.T ~ T
        r17.T ~ T
        r18.T ~ T
        r19.T ~ T
        r20.T ~ T
        r21.T ~ T
        r22.T ~ T
        r23.T ~ T
        r24.T ~ T
        r25.T ~ T
        r26.T ~ T
        r27.T ~ T
    end
end

@named eq = EquilibriumConstants()


"""
Create an equation system for this reaction based on the information in
Equation 5 and Table 2 of Fountoukis and Nenes (2007).

Rather than setting up for nonlinear rootfinding as is done in the original paper,
we create equations for a mass-action based solution where mass is moved between
the product and reactant species based on different between the current ratio and
the equilibrium ratio.

`activities` should be a system of equations for the activities of the species in this reaction
and `f_del` should be a dictionary of variables representing deliquescence fractions for each
species.
"""
function rxn_sys(r::Rxn, t, activities::ModelingToolkit.AbstractSystem, f_del)
    ap = speciesactivity(activities, r.product)
    ar = speciesactivity(activities, r.reactant)
    pv, ps = terms(r.product) # Product variables and stoichiometry coefficients
    rv, rs = terms(r.reactant) # Reactant variables and stoichiometry coefficients
    Dt = Differential(t)
    @constants T₀ = 293.15 [unit = u"K", description = "Standard temperature"]
    # These are the variables from Fountoukis and Nenes (2007) Table 2
    @constants K⁰ = r.K⁰ [
        unit = r.K⁰units,
        description = "Equilibrium constant at 298.15 K"
    ]
    @constants H_group = r.hgroup [description = "ΔH⁰ / (R * T₀) (unitless)"]
    @constants C_group = r.cgroup [description = "ΔC⁰ₚ / R (unitless)"]
    @variables K_eq(t) [unit = r.K⁰units, description = "Equilibrium constant"]
    @variables a_ratio(t) [unit = r.K⁰units, description = "Equilibrium constant"]
    @variables rawrate(t) [unit = r.K⁰units, description = "Pseudo reaction rate"]
    @variables rate2(t) [unit = r.K⁰units, description = "Normalized Pseudo reaction rate"]
    @variables rate(t) [unit = r.K⁰units, description = "Normalized Pseudo reaction rate"]
    @constants rateconst = 1e-9 [
        unit = r.K⁰units,
        description = "Rate constant (chosen to manage stiffness)"
    ]
    @constants zerorate = 0 [unit = r.K⁰units, description = "Zero rate"]
    @variables present(t) = 1 [
        description = "Whether the reactant is present (only used when reactant is a solid)",
    ]
    @constants unitconc = 1 [unit = u"mol/m_air^3", description = "Unit concentration"]
    @constants ratefactor = r.ratefactor [
        description = "Reaction-specific rate factor to manage stiffness",
    ]
    units = ([], [])
    for (i, vv) in enumerate((pv, rv))
        for v in vv
            x = Symbol(:unit, Symbolics.tosymbol(v, escape = false))
            push!(
                units[i],
                only(
                    @constants $x = 0.01 [
                    unit = ModelingToolkit.get_unit(v/rateconst),
                    description = "Unit conversion factor"
                ]
                )
            )
        end
    end
    #if (typeof(r.product) <: SaltLike) && (r.reactant isa Solid) # Deliquescence
    #    ar *= f_del[r.product]
    #end

    sys = ODESystem(
        [
         # Equation 5: Equilibrium constant
         K_eq ~ K⁰ * exp(-H_group * (T₀ / T - 1) - C_group * (1 + log(T₀ / T) - T₀ / T))
         a_ratio ~ ap / ar
         rawrate ~ (a_ratio - K_eq) * ratefactor
         # Decay reaction rate to zero as the reactant approaches zero concentration.
         rate2 ~ ifelse(
             rawrate > zerorate,
             min([rawrate; pv ./ units[1]]...),
             -min([-rawrate; rv ./ units[2]]...)
         )
         # Clip tiny reaction rates.
         rate ~ ifelse(abs(rate2) > rateconst, rate2, zerorate)],
        t,
        [K_eq, a_ratio, rawrate, rate2, rate],
        [];
        name = r.name
    )

    # Equations to move toward equilibrium
    ode_eqs = Dict()
    for (v, s, sign) in
        zip(vcat(pv, rv), vcat(ps, rs), vcat(fill(-1, length(pv)), fill(1, length(rv))))
        x = Symbol(r.name, :conv_, Symbolics.tosymbol(v, escape = false))
        conv = only(
            @constants $x = 1 [
            unit = ModelingToolkit.get_unit(v / rate / t),
            description = "Unit conversion factor"
        ]
        )
        ode_eqs[Dt(v)] = sign * sys.rate * s * conv
    end
    sys, ode_eqs
end

"""
Return the variable from the given system of activities that represents the activity of the given species.
"""
function speciesactivity(activities::ModelingToolkit.AbstractSystem, s::Species)
    vars = states(activities)
    varnames = [string(Symbolics.tosymbol(v, escape = false)) for v in vars]
    n = "a_$(nameof(s))"
    i = findfirst(isequal(n), varnames)
    if isnothing(i)
        error("Could not find activity variable for $n in $varnames")
    end
    ParentScope(vars[i])
end

"""
The activity of multiple species is the product of their activities
as shown in Table 2 of Fountoukis and Nenes (2007).
"""
function speciesactivity(activities::ModelingToolkit.AbstractSystem, s::AbstractVector)
    reduce(*, speciesactivity.((activities,), s))
end

"""
generate reactions, where slt, g, sld, and i are dictionaries of aqueous salts, gases, solids, and ions, respectively,
and H2Oaq is water.
"""
function generate_reactions(slt, g, sld, i, H2Oaq)
    [
     # NOTE: Assuming that H_group and C_group are zero when they are left out of Table 2.
     Rxn(
         sld[:CaNO32],
         slt[:CaNO32],
         6.067e5,
         u"mol^3/kg_water^3",
         -11.299,
         0.0,
         1e-12,
         :rxn1
     )
     Rxn(
         sld[:CaCl2],
         slt[:CaCl2],
         7.974e11,
         u"mol^3/kg_water^3",
         -14.087,
         0.0,
         1e-18,
         :rxn2
     )
     Rxn(sld[:CaSO4], slt[:CaSO4], 4.319e-5, u"mol^2/kg_water^2", 0.0, 0.0, 1e-2, :rxn3)
     Rxn(
         sld[:K2SO4],
         slt[:K2SO4],
         1.569e-2,
         u"mol^3/kg_water^3",
         -9.589,
         45.807,
         1e-10,
         :rxn4
     )
     Rxn(
         sld[:KHSO4],
         slt[:KHSO4],
         24.016,
         u"mol^2/kg_water^2",
         -8.423,
         17.964,
         1e-11,
         :rxn5
     )
     Rxn(sld[:KNO3], slt[:KNO3], 0.872, u"mol^2/kg_water^2", 14.075, 19.388, 1e-11, :rxn6)
     Rxn(sld[:KCl], slt[:KCl], 8.680, u"mol^2/kg_water^2", -6.167, 19.953, 1e-8, :rxn7)
     Rxn(sld[:MgSO4], slt[:MgSO4], 1.079e5, u"mol^2/kg_water^2", 36.798, 0.0, 1e-12, :rxn8)
     Rxn(
         sld[:MgNO32],
         slt[:MgNO32],
         2.507e15,
         u"mol^3/kg_water^3",
         -8.754,
         0.0,
         1e-22,
         :rxn9
     )
     Rxn(
         sld[:MgCl2],
         slt[:MgCl2],
         9.557e21,
         u"mol^3/kg_water^3",
         -1.347,
         0.0,
         1e-28,
         :rxn10
     )
     Rxn(i[:HSO4], [i[:H], i[:SO4]], 1.015e-2, u"mol/kg_water", 8.85, 25.14, 1e-5, :rxn11)
     Rxn(
         g[:NH3],
         i[:NH3],
         5.764e1,
         u"mol/kg_water/Constants.atm",
         13.79,
         -5.39,
         1e-12,
         :rxn12
     )
     Rxn(
         [i[:NH3], H2Oaq],
         [i[:NH4], i[:OH]],
         1.805e-5,
         u"mol/kg_water",
         -1.50,
         26.92,
         1e-6,
         :rxn13
     )
     Rxn(
         g[:HNO3],
         slt[:HNO3],
         2.511e6,
         u"mol^2/kg_water^2/Constants.atm",
         29.17,
         16.83,
         1e-13,
         :rxn14
     )
     Rxn(
         g[:HNO3],
         i[:HNO3],
         2.1e5,
         u"mol/kg_water/Constants.atm",
         29.17,
         16.83,
         1e-14,
         :rxn15
     )
     Rxn(
         g[:HCl],
         [i[:H], i[:Cl]],
         1.971e6,
         u"mol^2/kg_water^2/Constants.atm",
         30.20,
         19.91,
         1e-12,
         :rxn16
     )
     Rxn(
         g[:HCl],
         i[:HCl],
         2.5e3,
         u"mol/kg_water/Constants.atm",
         30.20,
         19.91,
         1e-14,
         :rxn17
     )
     Rxn(
         H2Oaq,
         [i[:H], i[:OH]],
         1.010e-14,
         u"mol^2/kg_water^2",
         -22.52,
         26.92,
         1e-8,
         :rxn18
     )
     Rxn(
         sld[:Na2SO4],
         slt[:Na2SO4],
         4.799e-1,
         u"mol^3/kg_water^3",
         0.98,
         39.75,
         1e-6,
         :rxn19
     )
     Rxn(
         sld[:NH42SO4],
         slt[:NH42SO4],
         1.87e0,
         u"mol^3/kg_water^3",
         -2.65,
         38.57,
         1e-7,
         :rxn20
     )
     Rxn(
         sld[:NH4Cl],
         [g[:NH3], g[:HCl]],
         1.086e-16,
         u"Constants.atm^2",
         -71.00,
         2.40,
         1e10,
         :rxn21
     )
     Rxn(sld[:NaNO3], slt[:NaNO3], 1.197e1,
         u"mol^2/kg_water^2", -8.22, 16.01, 1e-8, :rxn22)
     Rxn(sld[:NaCl], slt[:NaCl], 3.766e1, u"mol^2/kg_water^2", -1.56, 16.90, 1e-8, :rxn23)
     Rxn(
         sld[:NaHSO4],
         slt[:NaHSO4],
         2.413e4,
         u"mol^2/kg_water^2",
         0.79,
         14.75,
         1e-11,
         :rxn24
     )
     Rxn(
         sld[:NH4NO3],
         [g[:NH3], g[:HNO3]],
         4.199e-17,
         u"Constants.atm^2",
         -74.375,
         6.025,
         5e10,
         :rxn25
     )
     Rxn(
         sld[:NH4HSO4],
         slt[:NH4HSO4],
         1.383e0,
         u"mol^2/kg_water^2",
         -2.87,
         15.83,
         1e-7,
         :rxn26
     )
     Rxn(
         sld[:NH43HSO42],
         slt[:NH43HSO42],
         2.972e1,
         u"mol^5/kg_water^5",
         -5.19,
         54.40,
         1e-8,
         :rxn27
     )]
end

"""
Return a vector of the Species that are products or reactants in the given reactions.
"""
function active_species(rxns::Vector{Rxn})
    unique(vcat(vcat([[r.reactant, r.product] for r in rxns]...)...))
end
