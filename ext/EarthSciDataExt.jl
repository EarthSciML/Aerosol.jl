module EarthSciDataExt
using Aerosol, EarthSciData, ModelingToolkit, EarthSciMLBase, DynamicQuantities

@register_unit ppb 1

function EarthSciMLBase.couple2(
        c::Aerosol.ElementalCarbonCoupler,
        e::EarthSciData.NEI2016MonthlyEmisCoupler
    )
    c, e = c.sys, e.sys

    @constants(
        MW_C = 12.011e-3,
        [unit = u"kg/mol", description = "Carbon molar mass"],
        MW_air = 28.97e-3,
        [unit = u"kg/mol", description = "Dry air molar mass"],
        nmolpermol = 1.0e9,
        [unit = u"ppb", description = "nmol/mol, Conversion factor from mol to nmol"],
    )

    # In EarthSciData v0.15, emissions are in units of "kg/kg/s" (mass mixing ratio per second).
    # To convert to "ppb/s" (= nmol/mol/s):
    #   ppb/s = (kg_EC/kg_air/s) × (MW_air/MW_C) × nmolpermol
    # MW_air/MW_C converts mass mixing ratio to molar mixing ratio,
    # and nmolpermol (=1e9) converts mol/mol to nmol/mol (=ppb).
    uconv = nmolpermol * MW_air / MW_C
    return operator_compose(c, e, Dict(c.EC => e.PEC => uconv))
end

function EarthSciMLBase.couple2(
        c::Aerosol.ElementalCarbonCoupler,
        g::EarthSciData.GEOSFPCoupler
    )
    c, g = c.sys, g.sys
    c = param_to_var(c, :T, :P)
    return ConnectorSystem([c.T ~ g.I3₊T, c.P ~ g.P], c, g)
end

end
