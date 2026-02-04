module EarthSciDataExt
using Aerosol, EarthSciData, ModelingToolkit, EarthSciMLBase, DynamicQuantities

# Note: ppb unit is already registered in the main Aerosol module

function EarthSciMLBase.couple2(
        c::Aerosol.ElementalCarbonCoupler,
        e::EarthSciData.NEI2016MonthlyEmisCoupler
)
    c, e = c.sys, e.sys

    @constants(MW_C=12.011e-3,
        [unit=u"kg/mol", description="Carbon molar mass"],
        nmolpermol=1e9,
        [unit=u"ppb", description="nmol/mol, Conversion factor from mol to nmol"],
        MW_Air=28.97e-3,
        [unit=u"kg/mol", description="Molar mass of air"],)

    # Emissions are in units of "kg/kg air/s" and need to be converted to "ppb/s" or "nmol/mol/s".
    uconv = nmolpermol * MW_Air # Conversion factor with MW factored out.
    operator_compose(c, e, Dict(c.EC => e.PEC => uconv / MW_C))
end

function EarthSciMLBase.couple2(
        c::Aerosol.ElementalCarbonCoupler,
        g::EarthSciData.GEOSFPCoupler
)
    c, g = c.sys, g.sys
    c = param_to_var(c, :T, :P)
    ConnectorSystem([c.T ~ g.I3₊T, c.P ~ g.P], c, g)
end

end
