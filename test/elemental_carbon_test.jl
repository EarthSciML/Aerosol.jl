using Aerosol
using EarthSciData, EarthSciMLBase
using ModelingToolkit
using Dates
using Test

@testset "elemental carbon" begin
    domain = DomainInfo(DateTime(2016, 2, 3), DateTime(2016, 2, 4);
        latrange=deg2rad(40.0f0):deg2rad(2):deg2rad(44.0f0),
        lonrange=deg2rad(-97.0f0):deg2rad(2.5):deg2rad(-92.0f0),
        levrange=1:1, dtype=Float64)

    structural_simplify(ElementalCarbon())

    model = couple(
        NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain),
        GEOSFP("4x5", domain),
        ElementalCarbon()
    )

    sys = convert(ODESystem, model, prune=false)

    eqs = string(equations(sys))
    @test occursin("ElementalCarbon₊NEI2016MonthlyEmis_PEC(t)", eqs)

    obs = string(observed(sys))
    @test occursin("ElementalCarbon₊NEI2016MonthlyEmis_PEC(t) ~ (ElementalCarbon₊T*R*nmolpermol*NEI2016MonthlyEmis₊PEC(t)) / (ElementalCarbon₊P*MW_C)", obs)
    @test occursin("ElementalCarbon₊T(t) ~ GEOSFP₊I3₊T", obs)
end
