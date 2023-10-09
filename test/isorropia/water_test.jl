@test ModelingToolkit.get_unit(W_eq16) == u"kg_water/m_air^3"

@testset "water content" begin
    mw = Dict(Na_aq => 22.989769, SO4_aq => 96.0636, SO4_g => 96.0636, NH3_aq => 17.03052, NH3_g => 17.03052, NO3_aq => 62.0049, Cl_aq => 35.453,
        Ca_aq => 40.078, K_aq => 39.0983, Mg_aq => 24.305, H_aq => 1.00784, NH4_aq => 18.04,
        K2SO4_s => 174.259, KNO3_s => 101.1032, CaNO32_s => 164.1, HNO3_g => 63.01, HNO3_aq => 63.01) # g/mol

    ics = Dict([Na_aq => 0, SO4_aq => 10, NH3_aq => 3.4, NO3_aq => 2, Cl_aq => 0,
        Ca_aq => 0.4, K_aq => 0.33, Mg_aq => 0.0]) # ug/m3
    ics = Dict([k => ics[k] / 1e6 / mw[k] for k ∈ keys(ics)]) # ug/m3 / (1e6 ug/g) / g/mol = mol/m3
    ics[H_aq] = 2 * ics[SO4_aq] + ics[NO3_aq] + ics[Cl_aq]
    w2 = ModelingToolkit.subs_constants(W_eq16)
    w3 = ModelingToolkit.substitute(w2, ics)
    RHs = [10, 25, 40, 55, 65, 70, 75, 80, 85, 90] ./ 100.0
    ws = [Symbolics.value(ModelingToolkit.substitute(w3, [RH => v])) for v ∈ RHs] .* 1e9  # kg / m3 * (1e9 ug/kg) = ug/m3

    # TODO(CT): This should be larger than Fountoukis and Nenes (2007) Figure 6a once the rest of the salts are added in.
    ws_want = [9.274415859681406, 10.232660512301191, 11.34253892156881, 10.55022926762637, 12.790288404478733,
        13.626163770467416, 14.659016151077001, 16.101918073323315, 18.42839234123216, 23.167762687858897]
    @test ws ≈ ws_want
end

# Test that the activity of water is equal to the relative humidity.
test_subs = Dict([RH => 0.5])
@test ModelingToolkit.substitute(activity(H2Oaq), test_subs) == 0.5