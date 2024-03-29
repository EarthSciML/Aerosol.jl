using Aerosol
using NonlinearSolve, ModelingToolkit

@testset "VBS" begin
    Ci = [2.5, 1.8, 4.0, 4.0, 5.8, 4.8, 6.3, 8.0]
    ns = VBS(Ci)
    simplens = structural_simplify(ns)
    @unpack C_OA = simplens
    guess = [29]
    params = []
    prob = NonlinearProblem(simplens, guess, params)
    sol = solve(prob, TrustRegion())
    @test sol[C_OA] ≈ 10.6 atol=1e-2
end