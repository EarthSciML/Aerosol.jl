"""
A gas with the given partial pressure.
From Section 2.2 in Fountoukis and Nenes (2007), the activity of a gas is its partial pressure (in atm).
"""
@component function Gas(; name=:Gas, M_guess=nothing)
    @constants begin
        R = 8.31446261815324, [unit = u"m^3*Pa/K/mol", description = "Univ. gas constant"]
    end
    @constants begin
        p_one = 101325, [unit = u"Pa", description = "One atmosphere in Pascals"]
    end
    @variables begin
        T(t), [description = "Temperature", unit = u"K", guess = 293.15]
        p(t), [description = "Partial pressure", unit = u"Pa", guess = 1.0e-3]
        logp(t),
            [description = "Log of the partial pressure (in atm.) (dimensionless)", guess = log(1.0e-8)]
        M(t),
            [description = "Molarity of the gas in air", unit = u"mol/m^3", guess = M_guess]
    end
    eqs = [
        p ~ exp(logp) * p_one,
        M ~ p / R / T,
    ]
    return System(eqs, t; name)
end

"Gases in Isorropia II."
@component function Gases(; name=:Gases)
    @variables begin
        T(t), [unit = u"K", description = "Temperature", guess = 293.15]
    end
    HNO3 = Gas(; name=:HNO3, M_guess=6.0e-8)
    HCl = Gas(; name=:HCl, M_guess=1.6e-7)
    NH3 = Gas(; name=:NH3, M_guess=4.7e-8)
    eqs = [
        HNO3.T ~ T,
        HCl.T ~ T,
        NH3.T ~ T,
    ]
    return System(eqs, t; systems=[HNO3, HCl, NH3], name)
end
