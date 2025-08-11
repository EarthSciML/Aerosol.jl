@testitem "VBS_typical" begin
    using NonlinearSolve, ModelingToolkit
    using Aerosol
    Ci = [2.5, 1.8, 4.0, 4.0, 5.8, 4.8, 6.3, 8.0]
    ns = VBS(Ci)
    simplens = mtkcompile(ns)
    guess = [29]
    prob = NonlinearProblem(simplens, guess, [])
    sol = solve(prob, TrustRegion())
    @test sol[simplens.C_OA] ≈ 10.6 atol=1e-2
end

@testitem "different_temperature" begin
    using NonlinearSolve, ModelingToolkit
    Ci = [2.5, 1.8, 4.0, 4.0, 5.8, 4.8, 6.3, 8.0]
    ns = VBS(Ci)
    simplens = mtkcompile(ns)
    guess = [29]
    prob = NonlinearProblem(simplens, guess, [simplens.T=>285])
    sol = solve(prob, TrustRegion())
    @test sol[simplens.C_OA] ≈ 15.93 atol=1e-2
end
