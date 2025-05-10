module EarthSciDataExt
using Aerosol, EarthSciData, ModelingToolkit, EarthSciMLBase, DynamicQuantities

@register_unit ppb 1

function EarthSciMLBase.couple2(c::Aerosol.ElementalCarbonCoupler, e::EarthSciData.NEI2016MonthlyEmisCoupler)
    c, e = c.sys, e.sys

    @constants(
        MW_C = 12.011e-3, [unit = u"kg/mol", description="Carbon molar mass"],
        nmolpermol = 1e9, [unit = u"ppb", description="nmol/mol, Conversion factor from mol to nmol"],
        R = 8.31446261815324, [unit = u"m^3*Pa/mol/K", description="Ideal gas constant"],
    )

    # Emissions are in units of "kg/m3/s" and need to be converted to "ppb/s" or "nmol/mol/s".
    # To do this we need to convert kg of emissions to nmol of emissions,
    # and we need to convert m3 of air to mol of air.
    # nmol_emissions = kg_emissions * gperkg / MW_emission * nmolpermol = kg / kg/mol * nmol/mol = nmol
    # mol_air = m3_air / R / T * P = m3 / (m3*Pa/mol/K) / K * Pa = mol
    # So, the overall conversion is:
    # nmol_emissions / mol_air = (kg_emissions / MW_emission * nmolpermol) / (m3_air / R / T * P)
    uconv =  nmolpermol * R * c.T / c.P # Conversion factor with MW factored out.
    operator_compose(c, e, Dict(
        c.EC => e.PEC => uconv / MW_C,
    ))
end

function EarthSciMLBase.couple2(c::Aerosol.ElementalCarbonCoupler, g::EarthSciData.GEOSFPCoupler)
    c, g = c.sys, g.sys
    c = param_to_var(c, :T, :P)
    ConnectorSystem([
            c.T ~ g.I3₊T,
            c.P ~ g.P,
        ], c, g)
end

end
