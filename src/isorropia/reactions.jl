#@constants tinyrate = 1.e-99 [unit = u"mol/m_air^3/s", description = "Tiny reaction rate constant negative concentrations"]

"""
Define a reaction based on information in Table 2 of Fountoukis and Nenes (2007).

The left-hand side of the reaction equation is the equilibrium constant and the right-hand side is the ratio 
of the product and reactant activities.
"""
struct Rxn
    reactant
    product
    K⁰::Number
    K⁰units
    hgroup::Number
    cgroup::Number
    name
end

"""
Create an equation system for this reaction.
"""
function rxn_sys(r::Rxn, activities::ModelingToolkit.AbstractSystem)
    @constants T₀ = 293.15 [unit = u"K", description = "Standard temperature"]
    # These are the variables from Fountoukis and Nenes (2007) Table 2
    @constants K⁰ = r.K⁰ [unit = r.K⁰units, description = "Equilibrium constant at 298.15 K"]
    @constants H_group = r.hgroup [description = "ΔH⁰ / (R * T₀) (unitless)"]
    @constants C_group = r.cgroup [description = "ΔC⁰ₚ / R (unitless)"]
    @variables K_eq [unit = r.K⁰units, description = "Equilibrium constant"]
    ap = speciesactivity(activities, r.product)
    ar = speciesactivity(activities, r.reactant)
    @constants tinyap = 1e-20 [unit = ModelingToolkit.get_unit(ap), description = "Tiny activity to avoid negative concentration"]
    @constants tinyar = 1e-20 [unit = ModelingToolkit.get_unit(ar), description = "Tiny activity to avoid negative concentration"]
    sys = NonlinearSystem([
            # Equation 5: Equilibrium constant
            K_eq ~ K⁰ * exp(-H_group * (T₀ / T - 1) - C_group * (1 + log(T₀ / T) - T₀ / T))
            K_eq ~ max(ap, tinyap) / max(ar, tinyar)
        ], [K_eq], [K⁰, H_group, C_group, T₀, T₀₂, c_1, I_one, tinyap, tinyar]; name=r.name)
end


"""
Return the variable from the given system of activities that represents the activity of the given species.
"""
function speciesactivity(activities::ModelingToolkit.AbstractSystem, s::Species)
    vars = states(activities)
    varnames = [string(Symbolics.tosymbol(v, escape=false)) for v ∈ vars]
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
speciesactivity(activities::ModelingToolkit.AbstractSystem, s::AbstractVector) = reduce(*, speciesactivity.((activities,), s))



# NOTE: Assuming that H_group and C_group are zero when they are left out of Table 2.
reactions = [
    Rxn(CaNO32s, CaNO32_aqs, 6.067e5, u"mol^3/kg_water^3", -11.299, 0.0, :rxn1)
    Rxn(CaCl2s, CaCl2_aqs, 7.974e11, u"mol^3/kg_water^3", -14.087, 0.0, :rxn2)
    Rxn(CaSO4s, CaSO4_aqs, 4.319e-5, u"mol^2/kg_water^2", 0.0, 0.0, :rxn3)
    Rxn(K2SO4s, K2SO4_aqs, 1.569e-2, u"mol^3/kg_water^3", -9.589, 45.807, :rxn4)
    Rxn(KHSO4s, KHSO4_aqs, 24.016, u"mol^2/kg_water^2", -8.423, 17.964, :rxn5)
    Rxn(KNO3s, KNO3_aqs, 0.872, u"mol^2/kg_water^2", 14.075, 19.388, :rxn6)
    Rxn(KCls, KCl_aqs, 8.680, u"mol^2/kg_water^2", -6.167, 19.953, :rxn7)
    Rxn(MgSO4s, MgSO4_aqs, 1.079e5, u"mol^2/kg_water^2", 36.798, 0.0, :rxn8)
    Rxn(MgNO32s, MgNO32_aqs, 2.507e15, u"mol^3/kg_water^3", -8.754, 0.0, :rxn9)
    Rxn(MgCl2s, MgCl2_aqs, 9.557e21, u"mol^3/kg_water^3", -1.347, 0.0, :rxn10)
    Rxn(HSO4_ion, [H_ion, SO4_ion], 1.015e-2, u"mol/kg_water", 8.85, 25.14, :rxn11)
    Rxn(NH3g, NH3_ion, 5.764e1, u"mol/kg_water/atm", 13.79, -5.39, :rxn12)
    Rxn([NH3_ion, H2Oaq], [NH4_ion, OH_ion], 1.805e-5, u"mol/kg_water", -1.50, 26.92, :rxn13)
    Rxn(HNO3g, HNO3_aqs, 2.511e6, u"mol^2/kg_water^2/atm", 29.17, 16.83, :rxn14)
    Rxn(HNO3g, HNO3_ion, 2.1e5, u"mol/kg_water/atm", 29.17, 16.83, :rxn15)
    Rxn(HClg, [H_ion, Cl_ion], 1.971e6, u"mol^2/kg_water^2/atm", 30.20, 19.91, :rxn16)
    Rxn(HClg, HCl_ion, 2.5e3, u"mol/kg_water/atm", 30.20, 19.91, :rxn17)
    Rxn(H2Oaq, [H_ion, OH_ion], 1.010e-14, u"mol^2/kg_water^2", -22.52, 26.92, :rxn18)
    Rxn(Na2SO4s, Na2SO4_aqs, 4.799e-1, u"mol^3/kg_water^3", 0.98, 39.75, :rxn19)
    Rxn(NH42SO4s, NH42SO4_aqs, 1.87e0, u"mol^3/kg_water^3", -2.65, 38.57, :rxn20)
    Rxn(NH4Cls, [NH3g, HClg], 1.086e-16, u"atm^2", -71.00, 2.40, :rxn21)
    Rxn(NaNO3s, NaNO3_aqs, 1.197e1, u"mol^2/kg_water^2", -8.22, 16.01, :rxn22)
    Rxn(NaCls, NaCl_aqs, 3.766e1, u"mol^2/kg_water^2", -1.56, 16.90, :rxn23)
    Rxn(NaHSO4s, NaHSO4_aqs, 2.413e4, u"mol^2/kg_water^2", 0.79, 14.75, :rxn24)
    Rxn(NH4NO3s, [NH3g, HNO3g], 4.199e-17, u"atm^2", -74.375, 6.025, :rxn25)
    Rxn(NH4HSO4s, NH4HSO4_aqs, 1.383e0, u"mol^2/kg_water^2", -2.87, 15.83, :rxn26)
    Rxn(NH43HSO42s, NH43HSO42_aqs, 2.972e1, u"mol^5/kg_water^5", -5.19, 54.40, :rxn27)
]



"""
Return a vector of the Species that are products or reactants in the given reactions.
"""
active_species(rxns::Vector{Rxn}) = unique(vcat(vcat([[r.reactant, r.product] for r ∈ rxns]...)...))