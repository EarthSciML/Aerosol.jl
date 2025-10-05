@mtkmodel Gas begin
    @description """A gas with the given partial pressure.
    From Section 2.2 in Fountoukis and Nenes (2007), the activity of a gas is its partial pressure (in atm).
    """
    @constants begin
        R = 8.31446261815324, [unit = u"m^3*Pa/K/mol", description = "Univ. gas constant"]
        PaPerAtm = 101325, [unit = u"Pa/Constants.atm", description = "Unit conversion"]
    end
    @parameters begin
        T, [description = "Temperature", unit = u"K"]
    end
    @variables begin
        p(t), [description = "Partial pressure", unit = u"Constants.atm", guess=1e-4]
        M(t), [description = "Molarity of the gas in air", unit = u"mol/m^3", guess=1e-8]
    end
    @equations begin
        M ~ p * PaPerAtm / R / T
    end
end

@mtkmodel Gases begin
    @description "Gases in Isorropia II."
    @parameters begin
        T = 293.15, [unit = u"K", description = "Temperature"]
    end
    @components begin
        HNO3 = Gas()
        HCl = Gas()
        NH3 = Gas()
        H2SO4 = Gas()
    end
    @equations begin
        # From Section 3.3 (item 1) in Fountoukis and Nenes (2007), all H2SO4 immediately goes to aerosol phase.
        H2SO4.p ~ 0.0

        HNO3.T ~ T
        HCl.T ~ T
        NH3.T ~ T
        H2SO4.T ~ T
    end
end
