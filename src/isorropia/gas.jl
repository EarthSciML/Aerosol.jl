@mtkmodel Gas begin
    @description """A gas with the given partial pressure.
    From Section 2.2 in Fountoukis and Nenes (2007), the activity of a gas is its partial pressure (in atm).
    """
    @constants begin
        R = 8.31446261815324, [unit = u"m^3*Pa/K/mol", description = "Univ. gas constant"]
    end
    @parameters begin
        T, [description = "Temperature", unit = u"K"]
    end
    @structural_parameters begin
        M_guess=nothing
    end
    @constants begin
        p_one = 101325, [unit = u"Pa", description = "One atmosphere in Pascals"]
    end
    @variables begin
        p(t), [description = "Partial pressure", unit = u"Pa"]
        logp(t), [description = "Log of the partial pressure (in atm.)", guess=log(1e-8)]
        M(t), [description = "Molarity of the gas in air", unit = u"mol/m^3", guess=M_guess]
    end
    @equations begin
        p ~ exp(logp) * p_one
        M ~ p / R / T
    end
end

@mtkmodel Gases begin
    @description "Gases in Isorropia II."
    @parameters begin
        T = 293.15, [unit = u"K", description = "Temperature"]
    end
    @components begin
        HNO3 = Gas(M_guess=6e-8)
        HCl = Gas(M_guess=1.6e-7)
        NH3 = Gas(M_guess=4.7e-8)
    end
    @equations begin
        HNO3.T ~ T
        HCl.T ~ T
        NH3.T ~ T
    end
end
