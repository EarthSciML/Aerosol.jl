@testitem "elemental carbon" begin
    using Aerosol
    using EarthSciData, EarthSciMLBase
    using ModelingToolkit
    using Dates
    domain = DomainInfo(
        DateTime(2016, 2, 3),
        DateTime(2016, 2, 4);
        latrange = deg2rad(40.0f0):deg2rad(2):deg2rad(44.0f0),
        lonrange = deg2rad(-97.0f0):deg2rad(2.5):deg2rad(-92.0f0),
        levrange = 1:1
    )

    model = couple(
        NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain),
        GEOSFP("4x5", domain),
        ElementalCarbon()
    )

    sys = convert(System, model)

    eqs = string(equations(sys))
    @test occursin("ElementalCarbon₊NEI2016MonthlyEmis_PEC(t)", eqs)

    obs = string(observed(sys))
    @test occursin(
        "ElementalCarbon₊NEI2016MonthlyEmis_PEC(t) ~ (MW_air*nmolpermol*NEI2016MonthlyEmis₊PEC(t)) / MW_C",
        obs
    )
    @test occursin("ElementalCarbon₊T(t) ~ GEOSFP₊I3₊T(t)", obs)
end
