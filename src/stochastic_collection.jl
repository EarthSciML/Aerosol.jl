export StochasticCollectionCoalescence

using SymbolicUtils: @arrayop, @makearray

# Helper functions for @arrayop — take scalar symbolic arguments, return symbolic expressions.

function _sce_Z(n_dl, m_dl, xi_p)
    ifelse(n_dl > 0, xi_p * m_dl^2 / n_dl, 0.0)
end

function _sce_Q(n_dl, m_dl, xi_p)
    ifelse(n_dl > 0, xi_p^2 * m_dl^3 / n_dl^2, 0.0)
end

function _sce_R(n_dl, m_dl, xi_p)
    ifelse(n_dl > 0, xi_p^3 * m_dl^4 / n_dl^3, 0.0)
end

function _sce_f(n_dl, m_dl, xk)  # Eq. 13a + positivity Eq. 15
    xb = ifelse(n_dl > 0, m_dl / n_dl, xk)
    ifelse(xb < xk, 0.0,
        ifelse(xb > 2 * xk, 2 * n_dl / xk,
               4 * n_dl / xk - 2 * m_dl / xk^2))
end

function _sce_psi(n_dl, m_dl, xk)  # Eq. 13b + positivity Eq. 15
    xb = ifelse(n_dl > 0, m_dl / n_dl, xk)
    ifelse(xb < xk, 2 * n_dl / xk,
        ifelse(xb > 2 * xk, 0.0,
               2 * m_dl / xk^2 - 2 * n_dl / xk))
end

function _sce_alpha(n_dl, m_dl, xk)
    (_sce_f(n_dl, m_dl, xk) - _sce_psi(n_dl, m_dl, xk)) / (2 * xk)
end

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
SymbolicUtils.jl) for closure relations, linear distribution parameters, shifted arrays,
and summation reductions, with helper functions enabling complex logic within a single
@arrayop nesting level.

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

    # ---- Closure relations via @arrayop (level 1) ----
    Z_arr = @arrayop (k,) _sce_Z(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), XI_P) k in 1:I
    Q_arr = @arrayop (k,) _sce_Q(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), XI_P) k in 1:I

    # ---- Linear distribution parameters via @arrayop (level 1, Eq. 13a,b + Eq. 15a,b) ----
    f_arr = @arrayop (k,) _sce_f(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), x1 * 2.0^(k - 1)) k in 1:I
    ψ_arr = @arrayop (k,) _sce_psi(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), x1 * 2.0^(k - 1)) k in 1:I
    α_arr = @arrayop (k,) _sce_alpha(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), x1 * 2.0^(k - 1)) k in 1:I

    # ---- Shifted arrays via @makearray (level 1) ----
    N_prev_arr = @makearray N_prev_arr[1:I] begin
        N_prev_arr[1:1] => [0.0]
        N_prev_arr[2:I] => @arrayop (k,) Nk[k] * one_m3 k in 1:I-1
    end
    M_prev_arr = @makearray M_prev_arr[1:I] begin
        M_prev_arr[1:1] => [0.0]
        M_prev_arr[2:I] => @arrayop (k,) Mk[k] * (one_m3 / one_kg) k in 1:I-1
    end

    # ---- Tail sums via @arrayop (level 1) ----
    N_upper_arr = @arrayop (k,) ifelse(i >= k + 1, Nk[i] * one_m3, 0.0) k in 1:I i in 1:I
    M_lower_arr = @arrayop (k,) ifelse(i <= k - 1, Mk[i] * (one_m3 / one_kg), 0.0) k in 1:I i in 1:I

    # ---- Scalarize level-1 results for use in rate assembly ----
    N_dl_s = Symbolics.scalarize(@arrayop (k,) Nk[k] * one_m3 k in 1:I)
    M_dl_s = Symbolics.scalarize(@arrayop (k,) Mk[k] * (one_m3 / one_kg) k in 1:I)
    Z = Symbolics.scalarize(Z_arr)
    Q = Symbolics.scalarize(Q_arr)
    f = Symbolics.scalarize(f_arr)
    ψ = Symbolics.scalarize(ψ_arr)
    α = Symbolics.scalarize(α_arr)
    N_prev = Symbolics.scalarize(N_prev_arr)
    M_prev = Symbolics.scalarize(M_prev_arr)
    N_upper = Symbolics.scalarize(N_upper_arr)
    M_lower = Symbolics.scalarize(M_lower_arr)

    # ---- Kernel-specific rate computation ----
    if kernel_type == :constant
        dN_dl, dM_dl = _sce_rates_constant(N_dl_s, M_dl_s, Z, Q, f, ψ, α,
            N_prev, M_prev, N_upper, M_lower, x, I, K0_val)
    else  # :golovin
        R = Symbolics.scalarize(@arrayop (k,) _sce_R(Nk[k] * one_m3,
            Mk[k] * (one_m3 / one_kg), XI_P) k in 1:I)
        M_upper = Symbolics.scalarize(
            @arrayop (k,) ifelse(i >= k + 1, Mk[i] * (one_m3 / one_kg), 0.0) k in 1:I i in 1:I)
        Z_lower = Symbolics.scalarize(
            @arrayop (k,) ifelse(i <= k - 1,
                _sce_Z(Nk[i] * one_m3, Mk[i] * (one_m3 / one_kg), XI_P),
                0.0) k in 1:I i in 1:I)
        dN_dl, dM_dl = _sce_rates_golovin(N_dl_s, M_dl_s, Z, Q, R, f, ψ, α,
            N_prev, M_prev, N_upper, M_upper, M_lower, Z_lower, x, I, C_val)
    end

    # ---- Build equations with proper units ----
    eqs = vcat(
        [D(Nk[k]) ~ dN_dl[k] / (one_m3 * one_s) for k in 1:I],               # m⁻³ s⁻¹
        [D(Mk[k]) ~ dM_dl[k] * one_kg / (one_m3 * one_s) for k in 1:I],  # kg m⁻³ s⁻¹
    )

    return System(eqs, t; name, checks=ModelingToolkit.CheckComponents)
end

# ============================================================================
# RHS for constant kernel K(x,y) = K0
# Implements Eq. (9a,b) with K constant.
# ============================================================================
function _sce_rates_constant(N_dl, M_dl, Z, Q, f, ψ, α,
    N_prev, M_prev, N_upper, M_lower, x, I, K0)
    # S1+yS0 coefficients: a = 2xψ, b = f/2
    a_coeff = [2 * x[k] * ψ[k] for k in 1:I]
    b_coeff = [f[k] / 2 for k in 1:I]

    dN_dl = Vector{Any}(undef, I)
    dM_dl = Vector{Any}(undef, I)

    for k in 1:I
        # Shifted (k-1) quantities
        ψ_prev_k = k >= 2 ? ψ[k-1] : 0.0
        α_prev_k = k >= 2 ? α[k-1] : 0.0
        a_prev_k = k >= 2 ? a_coeff[k-1] : 0.0
        b_prev_k = k >= 2 ? b_coeff[k-1] : 0.0

        # ---- dN_k/dt (Eq. 9a) ----

        # Term 2: Σ_{i=1}^{k-2} S0(k-1 params, M[i], Z[i])
        term2_N = sum((ψ_prev_k * M_dl[i] + α_prev_k * Z[i]) for i in 1:max(k-2, 0); init=0.0)

        # Term 4: Σ_{i=1}^{k-1} S0(k params, M[i], Z[i])
        term4_N = sum((ψ[k] * M_dl[i] + α[k] * Z[i]) for i in 1:max(k-1, 0); init=0.0)

        # Assemble dN (N_upper and N_prev from @arrayop/@makearray)
        dN_dl[k] = K0 / 2 * N_prev[k]^2 + K0 * term2_N -
                   K0 * N_dl[k]^2 - K0 * term4_N -
                   K0 * N_dl[k] * N_upper[k]

        # ---- dM_k/dt (Eq. 9b) ----

        # Term 2: Σ_{i=1}^{k-2} (S1+yS0) using k-1 params
        term2_M = sum((a_prev_k * M_dl[i] + b_prev_k * Z[i] + α_prev_k * Q[i]) for i in 1:max(k-2, 0); init=0.0)

        # Term 4: Σ_{i=1}^{k-1} (S1+yS0) using k params
        term4_M = sum((a_coeff[k] * M_dl[i] + b_coeff[k] * Z[i] + α[k] * Q[i]) for i in 1:max(k-1, 0); init=0.0)

        # Assemble dM (M_lower, M_prev, N_upper from @arrayop/@makearray)
        dM_dl[k] = K0 * N_prev[k] * M_prev[k] + K0 * term2_M -
                   K0 * N_dl[k] * M_dl[k] - K0 * term4_M +
                   K0 * N_dl[k] * M_lower[k] - K0 * M_dl[k] * N_upper[k]
    end

    return dN_dl, dM_dl
end

# ============================================================================
# RHS for Golovin kernel K(x,y) = C*(x+y)
# Implements Eq. (9a,b) with K = C(x+y).
# ============================================================================
function _sce_rates_golovin(N_dl, M_dl, Z, Q, R, f, ψ, α,
    N_prev, M_prev, N_upper, M_upper, M_lower, Z_lower, x, I, C)
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
        Z_prev_k = k >= 2 ? Z[k-1] : 0.0
        a_prev_k = k >= 2 ? a_coeff[k-1] : 0.0
        b_prev_k = k >= 2 ? b_coeff[k-1] : 0.0
        α_prev_k = k >= 2 ? α[k-1] : 0.0
        d_prev_k = k >= 2 ? d_coeff[k-1] : 0.0
        e_prev_k = k >= 2 ? e_coeff[k-1] : 0.0
        g_prev_k = k >= 2 ? g_coeff[k-1] : 0.0

        # ---- dN_k/dt (Eq. 9a with K=C(x+y)) ----

        # Term 2: Σ_{i=1}^{k-2} (S1+yS0) using k-1 params
        term2_gN = sum((a_prev_k * M_dl[i] + b_prev_k * Z[i] + α_prev_k * Q[i]) for i in 1:max(k-2, 0); init=0.0)

        # Term 4: Σ_{i=1}^{k-1} (S1+yS0) using k params
        term4_gN = sum((a_coeff[k] * M_dl[i] + b_coeff[k] * Z[i] + α[k] * Q[i]) for i in 1:max(k-1, 0); init=0.0)

        # Assemble dN (N_prev, M_prev, N_upper, M_upper from @arrayop/@makearray)
        dN_dl[k] = C * N_prev[k] * M_prev[k] + C * term2_gN -
                   2 * C * N_dl[k] * M_dl[k] - C * term4_gN -
                   C * M_dl[k] * N_upper[k] - C * N_dl[k] * M_upper[k]

        # ---- dM_k/dt (Eq. 9b with K=C(x+y)) ----

        # Term 2: Σ_{i=1}^{k-2} (S2+yS1+y²S0) using k-1 params
        term2_gM = sum((d_prev_k * M_dl[i] + e_prev_k * Z[i] + g_prev_k * Q[i] + α_prev_k * R[i]) for i in 1:max(k-2, 0); init=0.0)

        # Term 4: Σ_{i=1}^{k-1} (S2+yS1+y²S0) using k params
        term4_gM = sum((d_coeff[k] * M_dl[i] + e_coeff[k] * Z[i] + g_coeff[k] * Q[i] + α[k] * R[i]) for i in 1:max(k-1, 0); init=0.0)

        # Assemble dM (N_prev, M_prev, Z_prev from @makearray/@arrayop; M_lower, Z_lower, N_upper, M_upper from @arrayop)
        dM_dl[k] = C * (N_prev[k] * Z_prev_k + M_prev[k]^2) + C * term2_gM -
                   C * (N_dl[k] * Z[k] + M_dl[k]^2) - C * term4_gM +
                   C * (M_dl[k] * M_lower[k] + N_dl[k] * Z_lower[k]) -
                   C * (Z[k] * N_upper[k] + M_dl[k] * M_upper[k])
    end

    return dN_dl, dM_dl
end
