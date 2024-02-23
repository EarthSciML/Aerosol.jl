ions = ISORROPIA.generate_ions(t)
salts = ISORROPIA.generate_salts(ions)
active_salts = collect(values(salts))
W = ISORROPIA.Water(t, active_salts)

@test ModelingToolkit.get_unit(equations(W)[1].rhs) == u"kg_water/m_air^3"

@testset "water content" begin
    mw = Dict(:Na => 22.989769, :SO4 => 96.0636, :NH3 => 17.03052, :NO3 => 62.0049, :Cl => 35.453,
        :Ca => 40.078, :K => 39.0983, :Mg => 24.305, :H => 1.00784, :NH4 => 18.04, :HSO4 => 97.064,
        :HNO3 => 63.01) # g/mol
    ics = Dict([:Na => 0, :SO4 => 10, :NH3 => 3.4, :NO3 => 2, :Cl => 0, :HSO4 => 0.0, :NH4 => 0.0,
        :Ca => 0.4, :K => 0.33, :Mg => 0.0]) # ug/m3
    ics = Dict([ions[k].m => ics[k] / 1e6 / mw[k] for k ∈ keys(ics)]) # ug/m3 / (1e6 ug/g) / g/mol = mol/m3
    ics[ions[:H].m] = 2 * ics[ions[:SO4].m] + ics[ions[:NO3].m] + ics[ions[:Cl].m]
    w2 = ModelingToolkit.subs_constants(equations(W)[1].rhs)
    w3 = ModelingToolkit.substitute(w2, ics)
    RHs = [10, 25, 40, 55, 65, 70, 75, 80, 85, 90] ./ 100.0
    ws = [Symbolics.value(ModelingToolkit.substitute(w3, [RH => v])) for v ∈ RHs] .* 1e9  # kg / m3 * (1e9 ug/kg) = ug/m3

    # TODO(CT): Is this correct?
    ws_want = [296.9078728183542, 299.3313075503013, 307.01470322554906, 317.03876229342495, 329.76309779878795, 337.1652639763639, 
        346.0988322675973, 358.08816326506627, 377.1477046922769, 415.7430650691068]
    @test ws ≈ ws_want
end

# Test that the activity of water is equal to the relative humidity.
#test_subs = Dict([RH => 0.5])
#@test ModelingToolkit.substitute(ISORROPIA.activity(ISORROPIA.H2O), test_subs) == 0.5