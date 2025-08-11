export ElementalCarbon

struct ElementalCarbonCoupler
    sys::Any
end

"""
A holder for elemental carbon concentrations. By default it does not include
any dynamics.

Default density is from: https://doi.org/10.5194/acp-23-4327-2023
"""
function ElementalCarbon(; name = :ElementalCarbon)
    params = @parameters begin
        ρ=1.1, [unit = u"kg/m^3", description = "Particle density"]
        d_p=0.8e-6,[unit = u"m", description = "Particle diameter"]
        T=298.15, [unit = u"K", description = "Temperature"]
        P=101325, [unit = u"Pa", description = "Pressure"]
    end
    vars = @variables begin
        EC(t)=1.0, [unit = u"ppb", description = "Elemental Carbon concentration"]
        v_placeholder(t),
        [
            unit = u"kg/m^2*K*Pa",
            description = "Placeholder to prevent parameters from being removed"
        ]
    end
    consts = @constants begin
        zero=0, [unit = ModelingToolkit.get_unit(EC / t), description = "Zero change"]
    end
    eqs = [
        D(EC) ~ zero,
        v_placeholder ~ ρ * d_p * T * P * 1e-20 # Meaningless equation to prevent parameters from being removed.
    ]
    System(
        eqs,
        t,
        vars,
       [params; consts];
        name = name,
        metadata = Dict(CoupleType => ElementalCarbonCoupler)
    )
end
