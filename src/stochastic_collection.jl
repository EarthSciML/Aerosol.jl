export StochasticCollectionCoalescence

using SymbolicUtils: @arrayop, @makearray

"""
    StochasticCollectionCoalescence(; name=:StochasticCollectionCoalescence, I=36, x1=1.6e-14, kernel_type=:constant, kernel_params=Dict(:K0 => 1e-10))

Solve the stochastic collection equation (SCE) for cloud droplet coalescence using
a two-moment method in discrete mass categories with `p = 2` (mass-doubling categories).

**Reference**: Tzivion, S., Feingold, G., and Levin, Z. (1987)
"An Efficient Numerical Solution to the Stochastic Collection Equation",
*Journal of the Atmospheric Sciences*, 44(21), 3139-3149.

The SCE describes the evolution of a drop size distribution due to gravitational
coalescence (collection). The continuous spectrum is divided into `I` discrete mass
categories with geometrically increasing boundaries (``x_{k+1} = 2 x_k``),
and two moments are tracked per category: number concentration ``N_k`` (m⁻³) and mass
concentration ``M_k`` (kg m⁻³). Closure is achieved via a nondimensional parameter ``\\bar{\\xi}_p``
that relates higher-order moments to these two (Eq. 8a,b).

The equations for each category follow Eq. (9a,b) in the paper, with incomplete
category integrals evaluated using a linear distribution function approximation
(Eq. 11-13) and positivity constraints (Eq. 15a,b).

The implementation uses symbolic array operations (`@arrayop` and `@makearray` from
SymbolicUtils.jl) to build transparent symbolic equations that the ModelingToolkit
system can fully analyze.

## Kernel Types

- `:constant` — ``K(x,y) = K_0``. Requires `kernel_params = Dict(:K0 => value)` in m³ s⁻¹.
- `:golovin` — ``K(x,y) = C(x+y)``. Requires `kernel_params = Dict(:C => value)` in m³ s⁻¹ kg⁻¹.

## Arguments

- `I`: Number of mass categories (default 36)
- `x1`: Lower mass boundary of the first category in kg (default ~1.6e-14 kg, corresponding to D≈3.125 μm)
- `kernel_type`: Collection kernel type, one of `:constant` or `:golovin`
- `kernel_params`: Dictionary of kernel parameters
"""
@component function StochasticCollectionCoalescence(;
    name=:StochasticCollectionCoalescence,
    I::Int=36, x1::Float64=1.6e-14,
    kernel_type::Symbol=:constant,
    kernel_params::Dict=Dict(:K0 => 1e-10))

    # ---- Category boundaries (Eq. 2, p=2) ----
    xb = [x1 * 2.0^(k - 1) for k in 1:(I+1)]

    # ---- Validate kernel ----
    if kernel_type == :constant
        K0_val = Float64(kernel_params[:K0])
    elseif kernel_type == :golovin
        C_val = Float64(kernel_params[:C])
    else
        error("Unsupported kernel_type: $kernel_type. Use :constant or :golovin.")
    end

    XI_P = 1.0625  # Closure parameter (Eq. 7, B10): 1 ≤ ξ_p ≤ 9/8, mean ≈ 1.0625

    @constants begin
        one_m3 = 1.0, [description = "Unit volume for dimensional analysis", unit = u"m^3"]
        one_kg = 1.0, [description = "Unit mass for dimensional analysis", unit = u"kg"]
        one_s = 1.0, [description = "Unit time for dimensional analysis", unit = u"s"]
    end

    @variables begin
        Nk(t)[1:I], [description = "Number concentration in category k", unit = u"m^-3"]
        Mk(t)[1:I], [description = "Mass concentration in category k", unit = u"kg*m^-3"]
    end

    x = xb[1:I]  # category lower boundaries

    # ---- Dimensionless quantities ----
    N_dl = Symbolics.scalarize(Nk .* one_m3)              # dimensionless
    M_dl = Symbolics.scalarize(Mk .* (one_m3 / one_kg))   # dimensionless

    # ---- Closure relations (Eq. 8a,b) ----
    Z = [ifelse(N_dl[k] > 0, XI_P * M_dl[k]^2 / N_dl[k], 0.0) for k in 1:I]
    Q = [ifelse(N_dl[k] > 0, XI_P^2 * M_dl[k]^3 / N_dl[k]^2, 0.0) for k in 1:I]

    # ---- Linear parameters (Eq. 13a,b with positivity Eq. 15a,b) ----
    x_bar = [ifelse(N_dl[k] > 0, M_dl[k] / N_dl[k], x[k]) for k in 1:I]
    f = [ifelse(x_bar[k] < x[k], 0.0,
          ifelse(x_bar[k] > 2 * x[k], 2 * N_dl[k] / x[k],
                 4 * N_dl[k] / x[k] - 2 * M_dl[k] / x[k]^2)) for k in 1:I]
    ψ = [ifelse(x_bar[k] < x[k], 2 * N_dl[k] / x[k],
          ifelse(x_bar[k] > 2 * x[k], 0.0,
                 2 * M_dl[k] / x[k]^2 - 2 * N_dl[k] / x[k])) for k in 1:I]

    # ---- S0 coefficient: α = (f-ψ)/(2x) ----
    α = [(f[k] - ψ[k]) / (2 * x[k]) for k in 1:I]

    # ---- Kernel-specific rate computation ----
    if kernel_type == :constant
        dN_dl, dM_dl = _sce_symbolic_constant(N_dl, M_dl, Z, Q, f, ψ, α, x, I, K0_val)
    else  # :golovin
        R = [ifelse(N_dl[k] > 0, XI_P^3 * M_dl[k]^4 / N_dl[k]^3, 0.0) for k in 1:I]
        dN_dl, dM_dl = _sce_symbolic_golovin(N_dl, M_dl, Z, Q, R, f, ψ, α, x, I, C_val)
    end

    # ---- Build equations with proper units ----
    eqs = vcat(
        [D(Nk[k]) ~ dN_dl[k] / (one_m3 * one_s) for k in 1:I],               # m⁻³ s⁻¹
        [D(Mk[k]) ~ dM_dl[k] * one_kg / (one_m3 * one_s) for k in 1:I],  # kg m⁻³ s⁻¹
    )

    return System(eqs, t; name, checks=ModelingToolkit.CheckComponents)
end

# ============================================================================
# Symbolic RHS for constant kernel K(x,y) = K0
# Implements Eq. (9a,b) with K constant, using symbolic array operations.
# ============================================================================
function _sce_symbolic_constant(N_dl, M_dl, Z, Q, f, ψ, α, x, I, K0)
    # S1+yS0 coefficients: a = 2xψ, b = f/2
    a_coeff = [2 * x[k] * ψ[k] for k in 1:I]
    b_coeff = [f[k] / 2 for k in 1:I]

    dN_dl = Vector{Any}(undef, I)
    dM_dl = Vector{Any}(undef, I)

    for k in 1:I
        # Shifted (k-1) quantities
        N_prev_k = k >= 2 ? N_dl[k-1] : 0.0
        M_prev_k = k >= 2 ? M_dl[k-1] : 0.0
        ψ_prev_k = k >= 2 ? ψ[k-1] : 0.0
        α_prev_k = k >= 2 ? α[k-1] : 0.0
        a_prev_k = k >= 2 ? a_coeff[k-1] : 0.0
        b_prev_k = k >= 2 ? b_coeff[k-1] : 0.0

        # ---- dN_k/dt (Eq. 9a) ----

        # Term 2: Cross-coag gain, Σ_{i=1}^{k-2} S0(f[k-1],ψ[k-1],x[k-1], M[i], Z[i])
        term2_N = sum((ψ_prev_k * M_dl[i] + α_prev_k * Z[i]) for i in 1:max(k-2, 0); init=0.0)

        # Term 4: Cross-coag loss, Σ_{i=1}^{k-1} S0(f[k],ψ[k],x[k], M[i], Z[i])
        term4_N = sum((ψ[k] * M_dl[i] + α[k] * Z[i]) for i in 1:max(k-1, 0); init=0.0)

        # Term 6: Upper tail sum, Σ_{i=k+1}^{I} N[i]
        N_upper = sum(N_dl[i] for i in (k+1):I; init=0.0)

        # Assemble dN: Terms 3+5 combined as -K0*N[k]^2
        dN_dl[k] = K0 / 2 * N_prev_k^2 + K0 * term2_N -
                   K0 * N_dl[k]^2 - K0 * term4_N -
                   K0 * N_dl[k] * N_upper

        # ---- dM_k/dt (Eq. 9b) ----

        # Term 2: Cross-coag mass gain, Σ_{i=1}^{k-2} (S1+yS0) using k-1 params
        term2_M = sum((a_prev_k * M_dl[i] + b_prev_k * Z[i] + α_prev_k * Q[i]) for i in 1:max(k-2, 0); init=0.0)

        # Term 4: Cross-coag mass loss, Σ_{i=1}^{k-1} (S1+yS0) using k params
        term4_M = sum((a_coeff[k] * M_dl[i] + b_coeff[k] * Z[i] + α[k] * Q[i]) for i in 1:max(k-1, 0); init=0.0)

        # Term 5: Lower cumulative sum of M, Σ_{i=1}^{k-1} M[i]
        M_lower = sum(M_dl[i] for i in 1:max(k-1, 0); init=0.0)

        # Assemble dM
        dM_dl[k] = K0 * N_prev_k * M_prev_k + K0 * term2_M -
                   K0 * N_dl[k] * M_dl[k] - K0 * term4_M +
                   K0 * N_dl[k] * M_lower - K0 * M_dl[k] * N_upper
    end

    return dN_dl, dM_dl
end

# ============================================================================
# Symbolic RHS for Golovin kernel K(x,y) = C*(x+y)
# Implements Eq. (9a,b) with K = C(x+y), using symbolic array operations.
# ============================================================================
function _sce_symbolic_golovin(N_dl, M_dl, Z, Q, R, f, ψ, α, x, I, C)
    # Coefficients
    a_coeff = [2 * x[k] * ψ[k] for k in 1:I]
    b_coeff = [f[k] / 2 for k in 1:I]
    d_coeff = [4 * x[k]^2 * ψ[k] for k in 1:I]
    e_coeff = [x[k] * f[k] / 2 + 2 * x[k] * ψ[k] for k in 1:I]
    g_coeff = [f[k] - ψ[k] for k in 1:I]

    dN_dl = Vector{Any}(undef, I)
    dM_dl = Vector{Any}(undef, I)

    for k in 1:I
        # Shifted (k-1) quantities
        N_prev_k = k >= 2 ? N_dl[k-1] : 0.0
        M_prev_k = k >= 2 ? M_dl[k-1] : 0.0
        Z_prev_k = k >= 2 ? Z[k-1] : 0.0
        a_prev_k = k >= 2 ? a_coeff[k-1] : 0.0
        b_prev_k = k >= 2 ? b_coeff[k-1] : 0.0
        α_prev_k = k >= 2 ? α[k-1] : 0.0
        d_prev_k = k >= 2 ? d_coeff[k-1] : 0.0
        e_prev_k = k >= 2 ? e_coeff[k-1] : 0.0
        g_prev_k = k >= 2 ? g_coeff[k-1] : 0.0

        # ---- dN_k/dt (Eq. 9a with K=C(x+y)) ----

        # Term 2: cross-coag gain, Σ_{i=1}^{k-2} (S1+yS0) using k-1 params
        term2_gN = sum((a_prev_k * M_dl[i] + b_prev_k * Z[i] + α_prev_k * Q[i]) for i in 1:max(k-2, 0); init=0.0)

        # Term 4: cross-coag loss, Σ_{i=1}^{k-1} (S1+yS0) using k params
        term4_gN = sum((a_coeff[k] * M_dl[i] + b_coeff[k] * Z[i] + α[k] * Q[i]) for i in 1:max(k-1, 0); init=0.0)

        # Upper tail sums
        N_upper = sum(N_dl[i] for i in (k+1):I; init=0.0)
        M_upper = sum(M_dl[i] for i in (k+1):I; init=0.0)

        # Assemble dN
        dN_dl[k] = C * N_prev_k * M_prev_k + C * term2_gN -
                   2 * C * N_dl[k] * M_dl[k] - C * term4_gN -
                   C * M_dl[k] * N_upper - C * N_dl[k] * M_upper

        # ---- dM_k/dt (Eq. 9b with K=C(x+y)) ----

        # Term 2: cross-coag mass gain, Σ_{i=1}^{k-2} (S2+yS1+y²S0) using k-1 params
        term2_gM = sum((d_prev_k * M_dl[i] + e_prev_k * Z[i] + g_prev_k * Q[i] + α_prev_k * R[i]) for i in 1:max(k-2, 0); init=0.0)

        # Term 4: cross-coag mass loss, Σ_{i=1}^{k-1} (S2+yS1+y²S0) using k params
        term4_gM = sum((d_coeff[k] * M_dl[i] + e_coeff[k] * Z[i] + g_coeff[k] * Q[i] + α[k] * R[i]) for i in 1:max(k-1, 0); init=0.0)

        # Lower sums
        M_lower = sum(M_dl[i] for i in 1:max(k-1, 0); init=0.0)
        Z_lower = sum(Z[i] for i in 1:max(k-1, 0); init=0.0)

        # Assemble dM
        dM_dl[k] = C * (N_prev_k * Z_prev_k + M_prev_k^2) + C * term2_gM -
                   C * (N_dl[k] * Z[k] + M_dl[k]^2) - C * term4_gM +
                   C * (M_dl[k] * M_lower + N_dl[k] * Z_lower) -
                   C * (Z[k] * N_upper + M_dl[k] * M_upper)
    end

    return dN_dl, dM_dl
end
