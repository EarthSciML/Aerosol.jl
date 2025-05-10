export ElementalCarbon

struct ElementalCarbonCoupler
    sys::Any
end

"""
A holder for elemental carbon concentrations. By default it does not include
any dynamics.

Default density is from: https://doi.org/10.5194/acp-23-4327-2023
"""
function ElementalCarbon(; name=:ElementalCarbon)
    params = @parameters(
        ρ = 1.1, [unit = u"kg/m^3", description = "Particle density"],
        d_p = 0.8e-6, [unit = u"m", description = "Particle diameter"],
        T = 298.15, [unit = u"K", description = "Temperature"],
        P = 101325, [unit = u"Pa", description = "Pressure"],
    )
    vars = @variables(
        EC(t), [unit = u"ppb", description = "Elemental Carbon concentration"],
        v_placeholder(t), [unit = u"kg/m^2*K*Pa", description = "Placeholder to prevent parameters from being removed"],
    )
    @constants zero = 0 [unit = ModelingToolkit.get_unit(EC/t), description = "Zero change"]
    eqs = [
        D(EC) ~ zero,
        v_placeholder ~ ρ * d_p * T * P * 1e-20, # Meaningless equation to prevent parameters from being removed.
    ]
    ODESystem(eqs, t, vars, params; name=name,
        metadata=Dict(:coupletype => ElementalCarbonCoupler))
end
