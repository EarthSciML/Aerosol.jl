export StochasticCollectionCoalescence

using Einsum
using ModelingToolkit.Symbolics: @register_array_symbolic

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

The implementation uses registered array symbolic functions to keep symbolic
expressions compact and avoid memory blowup for large numbers of categories.

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
        kernel_flag = 1.0  # constant kernel
        kernel_coeff = K0_val
    elseif kernel_type == :golovin
        C_val = Float64(kernel_params[:C])
        kernel_flag = 2.0  # Golovin kernel
        kernel_coeff = C_val
    else
        error("Unsupported kernel_type: $kernel_type. Use :constant or :golovin.")
    end

    @constants begin
        one_m3 = 1.0, [description = "Unit volume for dimensional analysis", unit = u"m^3"]
        one_kg = 1.0, [description = "Unit mass for dimensional analysis", unit = u"kg"]
        one_s = 1.0, [description = "Unit time for dimensional analysis", unit = u"s"]
    end

    @variables begin
        Nk(t)[1:I], [description = "Number concentration in category k", unit = u"m^-3"]
        Mk(t)[1:I], [description = "Mass concentration in category k", unit = u"kg*m^-3"]
    end

    # Build dimensionless state vector for the registered function.
    # The registered function operates on dimensionless quantities to avoid
    # unit-related issues in the symbolic system.
    N_dimless = Symbolics.scalarize(Nk .* one_m3)              # dimensionless
    M_dimless = Symbolics.scalarize(Mk .* (one_m3 / one_kg))   # dimensionless
    state = vcat(N_dimless, M_dimless)

    # Call the registered RHS function (returns dimensionless rates in s⁻¹ equivalents)
    result = _sce_rhs(state, xb, kernel_flag, kernel_coeff)

    # Build equations with proper units
    eqs = vcat(
        map(k -> D(Nk[k]) ~ result[k] / (one_m3 * one_s), 1:I),               # m⁻³ s⁻¹
        map(k -> D(Mk[k]) ~ result[I + k] * one_kg / (one_m3 * one_s), 1:I),  # kg m⁻³ s⁻¹
    )

    # Use CheckComponents to skip unit checking for registered array functions.
    # The unit checker cannot trace through opaque registered functions, but the
    # dimensional analysis is correct by construction (tested numerically).
    return System(eqs, t; name, checks=ModelingToolkit.CheckComponents)
end

# ============================================================================
# Registered array function for SCE right-hand side.
# Takes a state vector [N₁...Nᵢ, M₁...Mᵢ] (all dimensionless) plus
# category boundaries and kernel parameters.
# Returns [dN₁/dt...dNᵢ/dt, dM₁/dt...dMᵢ/dt] (dimensionless rates).
# ============================================================================
function _sce_rhs(state::AbstractVector, xb::Vector{Float64},
                  kernel_flag::Float64, kernel_coeff::Float64)
    n = length(state) ÷ 2
    N = collect(state[1:n])
    M = collect(state[n+1:2n])
    dN = zeros(eltype(state), n)
    dM = zeros(eltype(state), n)

    if kernel_flag == 1.0
        _sce_constant_rhs!(dN, dM, N, M, n, xb, kernel_coeff)
    elseif kernel_flag == 2.0
        _sce_golovin_rhs!(dN, dM, N, M, n, xb, kernel_coeff)
    end

    return vcat(dN, dM)
end

@register_array_symbolic _sce_rhs(state::AbstractVector, xb::Vector{Float64},
                                  kernel_flag::Float64, kernel_coeff::Float64) begin
    size = (length(state),)
    eltype = Real
    ndims = 1
end

# ============================================================================
# Closure parameter (Eq. 7, B10)
# For p=2: 1 ≤ ξ_p ≤ 9/8. Mean value ≈ 1.0625.
# ============================================================================
const _SCE_XI_P = 1.0625

# ============================================================================
# Linear approximation parameters (Eq. 13a,b with positivity Eq. 15a,b)
# Returns (fk, ψk) for the linear distribution function approximation
# in category k with lower boundary xk.
# ============================================================================
function _sce_linear_params(Nk_val, Mk_val, xk)
    # x̄_k = Mk/Nk (mean mass in category k)
    x_bar = Nk_val > 0 ? Mk_val / Nk_val : xk

    if x_bar < xk  # Eq. 15b: x̄_k < x_k → set fk=0 for positivity
        fk = 0.0
        ψk = 2 * Nk_val / xk
    elseif x_bar > 2 * xk  # Eq. 15a: x̄_k > x_{k+1} → set ψk=0 for positivity
        fk = 2 * Nk_val / xk
        ψk = 0.0
    else  # Normal case: Eq. 13a,b expanded
        # fk = (2Nk/xk)(2 - x̄k/xk) = 4Nk/xk - 2Mk/xk²
        fk = 4 * Nk_val / xk - 2 * Mk_val / xk^2
        # ψk = (2Nk/xk)(x̄k/xk - 1) = 2Mk/xk² - 2Nk/xk
        ψk = 2 * Mk_val / xk^2 - 2 * Nk_val / xk
    end
    return fk, ψk
end

# ============================================================================
# Closure relations (Eq. 8a,b)
# Z = ξ_p * M² / N  (second moment proxy)
# Q = ξ_p² * M³ / N²  (third moment proxy)
# R = ξ_p³ * M⁴ / N³  (fourth moment proxy, needed for Golovin mass eqs)
# ============================================================================
@inline _sce_Z(Mk, Nk) = Nk > 0 ? _SCE_XI_P * Mk^2 / Nk : 0.0
@inline _sce_Q(Mk, Nk) = Nk > 0 ? _SCE_XI_P^2 * Mk^3 / Nk^2 : 0.0
@inline _sce_R(Mk, Nk) = Nk > 0 ? _SCE_XI_P^3 * Mk^4 / Nk^3 : 0.0

# ============================================================================
# Incomplete integrals (Eq. 12) for linear approximation
# S0 = ∫_{x_{k+1}-y}^{x_{k+1}} n_m(x) dx  (number over incomplete interval)
# S1 = ∫ x·n_m(x) dx  (mass over incomplete interval)
# yS0 = ∫ y·n_m(x) dx  (y-weighted number over incomplete interval)
# ============================================================================
@inline _sce_S0(fm, ψm, xm, Mi, Zi) = ψm * Mi + (fm - ψm) / (2 * xm) * Zi
@inline _sce_S1(fm, ψm, xm, Mi, Zi) = 2 * xm * ψm * Mi + (fm - 2 * ψm) / 2 * Zi
@inline _sce_yS0(fm, ψm, xm, Zi, Qi) = ψm * Zi + (fm - ψm) / (2 * xm) * Qi

# ============================================================================
# Numerical RHS for constant kernel K(x,y) = K0
# Implements Eq. (9a,b) with K constant, using einsum array operations.
# ============================================================================
function _sce_constant_rhs!(dN, dM, N, M, I, xb, K0)
    T = eltype(N)
    x = xb[1:I]  # category lower boundaries

    # Precompute closure relations (Eq. 8a,b)
    @einsum Z[k] := _sce_Z(M[k], N[k])
    @einsum Q[k] := _sce_Q(M[k], N[k])

    # Linear params (Eq. 13/15) - broadcasting for tuple return
    fψ = _sce_linear_params.(N, M, x)
    f = getindex.(fψ, 1)
    ψ = getindex.(fψ, 2)

    # S0 coefficients: S0(f,ψ,x,Mi,Zi) = ψ*Mi + α*Zi where α = (f-ψ)/(2x)
    @einsum α[k] := (f[k] - ψ[k]) / (2 * x[k])

    # S1+yS0 coefficients: a*Mi + b*Zi + c*Qi where a=2xψ, b=f/2, c=α
    @einsum a_coeff[k] := 2 * x[k] * ψ[k]
    @einsum b_coeff[k] := f[k] / 2

    # Shifted (k-1) quantities for Term 2, which uses params from category k-1
    N_prev = vcat([zero(T)], N[1:I-1])
    M_prev = vcat([zero(T)], M[1:I-1])
    ψ_prev = vcat([zero(T)], ψ[1:I-1])
    α_prev = vcat([zero(T)], α[1:I-1])
    a_prev = vcat([zero(T)], a_coeff[1:I-1])
    b_prev = vcat([zero(T)], b_coeff[1:I-1])

    # Triangular mask matrices for index-bounded sums
    mask_lt2 = [i <= k - 2 ? one(T) : zero(T) for k in 1:I, i in 1:I]  # i=1..k-2
    mask_lt1 = [i <= k - 1 ? one(T) : zero(T) for k in 1:I, i in 1:I]  # i=1..k-1
    mask_gt  = [i >= k + 1 ? one(T) : zero(T) for k in 1:I, i in 1:I]  # i=k+1..I

    # ---- dN_k/dt (Eq. 9a) ----

    # Term 2: Cross-coag gain, Σ_{i=1}^{k-2} S0(f[k-1],ψ[k-1],x[k-1], M[i], Z[i])
    @einsum term2_N[k] := mask_lt2[k, i] * (ψ_prev[k] * M[i] + α_prev[k] * Z[i])

    # Term 4: Cross-coag loss, Σ_{i=1}^{k-1} S0(f[k],ψ[k],x[k], M[i], Z[i])
    @einsum term4_N[k] := mask_lt1[k, i] * (ψ[k] * M[i] + α[k] * Z[i])

    # Term 6: Upper tail sum, Σ_{i=k+1}^{I} N[i]
    @einsum N_upper[k] := mask_gt[k, i] * N[i]

    # Assemble dN: Terms 3+5 combined as -K0*N[k]^2
    @einsum dN[k] = (K0 / 2) * N_prev[k]^2 + K0 * term2_N[k] - K0 * N[k]^2 - K0 * term4_N[k] - K0 * N[k] * N_upper[k]

    # ---- dM_k/dt (Eq. 9b) ----

    # Term 2: Cross-coag mass gain, Σ_{i=1}^{k-2} (S1+yS0) using k-1 params
    @einsum term2_M[k] := mask_lt2[k, i] * (a_prev[k] * M[i] + b_prev[k] * Z[i] + α_prev[k] * Q[i])

    # Term 4: Cross-coag mass loss, Σ_{i=1}^{k-1} (S1+yS0) using k params
    @einsum term4_M[k] := mask_lt1[k, i] * (a_coeff[k] * M[i] + b_coeff[k] * Z[i] + α[k] * Q[i])

    # Term 5: Lower cumulative sum of M, Σ_{i=1}^{k-1} M[i]
    @einsum M_lower[k] := mask_lt1[k, i] * M[i]

    # Assemble dM
    @einsum dM[k] = K0 * N_prev[k] * M_prev[k] + K0 * term2_M[k] - K0 * N[k] * M[k] - K0 * term4_M[k] + K0 * N[k] * M_lower[k] - K0 * M[k] * N_upper[k]

    return nothing
end

# ============================================================================
# Numerical RHS for Golovin kernel K(x,y) = C*(x+y)
# Implements Eq. (9a,b) with K = C(x+y), using einsum array operations.
# ============================================================================
function _sce_golovin_rhs!(dN, dM, N, M, I, xb, C)
    T = eltype(N)
    x = xb[1:I]  # category lower boundaries

    # Precompute closure relations (Eq. 8a,b)
    @einsum Z[k] := _sce_Z(M[k], N[k])
    @einsum Q[k] := _sce_Q(M[k], N[k])
    @einsum R[k] := _sce_R(M[k], N[k])

    # Linear params (Eq. 13/15) - broadcasting for tuple return
    fψ = _sce_linear_params.(N, M, x)
    f = getindex.(fψ, 1)
    ψ = getindex.(fψ, 2)

    # S0 coefficient: α = (f-ψ)/(2x)
    @einsum α[k] := (f[k] - ψ[k]) / (2 * x[k])

    # S1+yS0 coefficients (for dN): a=2xψ, b=f/2, c=α
    @einsum a_coeff[k] := 2 * x[k] * ψ[k]
    @einsum b_coeff[k] := f[k] / 2

    # S2+yS1+y2S0 coefficients (for dM): d=4x²ψ, e=xf/2+2xψ, g=f-ψ, h=α
    @einsum d_coeff[k] := 4 * x[k]^2 * ψ[k]
    @einsum e_coeff[k] := x[k] * f[k] / 2 + 2 * x[k] * ψ[k]
    @einsum g_coeff[k] := f[k] - ψ[k]

    # Shifted (k-1) quantities for Term 2
    N_prev = vcat([zero(T)], N[1:I-1])
    M_prev = vcat([zero(T)], M[1:I-1])
    Z_prev = vcat([zero(T)], Z[1:I-1])
    a_prev = vcat([zero(T)], a_coeff[1:I-1])
    b_prev = vcat([zero(T)], b_coeff[1:I-1])
    α_prev = vcat([zero(T)], α[1:I-1])
    d_prev = vcat([zero(T)], d_coeff[1:I-1])
    e_prev = vcat([zero(T)], e_coeff[1:I-1])
    g_prev = vcat([zero(T)], g_coeff[1:I-1])

    # Triangular mask matrices for index-bounded sums
    mask_lt2 = [i <= k - 2 ? one(T) : zero(T) for k in 1:I, i in 1:I]  # i=1..k-2
    mask_lt1 = [i <= k - 1 ? one(T) : zero(T) for k in 1:I, i in 1:I]  # i=1..k-1
    mask_gt  = [i >= k + 1 ? one(T) : zero(T) for k in 1:I, i in 1:I]  # i=k+1..I

    # ---- dN_k/dt (Eq. 9a with K=C(x+y)) ----

    # Term 2: cross-coag gain, Σ_{i=1}^{k-2} (S1+yS0) using k-1 params
    @einsum term2_gN[k] := mask_lt2[k, i] * (a_prev[k] * M[i] + b_prev[k] * Z[i] + α_prev[k] * Q[i])

    # Term 4: cross-coag loss, Σ_{i=1}^{k-1} (S1+yS0) using k params
    @einsum term4_gN[k] := mask_lt1[k, i] * (a_coeff[k] * M[i] + b_coeff[k] * Z[i] + α[k] * Q[i])

    # Upper tail sums for Term 6
    @einsum N_upper[k] := mask_gt[k, i] * N[i]
    @einsum M_upper[k] := mask_gt[k, i] * M[i]

    # Assemble dN: Terms 3+5 combined as -2C*N[k]*M[k]
    @einsum dN[k] = C * N_prev[k] * M_prev[k] + C * term2_gN[k] - 2 * C * N[k] * M[k] - C * term4_gN[k] - C * M[k] * N_upper[k] - C * N[k] * M_upper[k]

    # ---- dM_k/dt (Eq. 9b with K=C(x+y)) ----

    # Term 2: cross-coag mass gain, Σ_{i=1}^{k-2} (S2+yS1+y2S0) using k-1 params
    @einsum term2_gM[k] := mask_lt2[k, i] * (d_prev[k] * M[i] + e_prev[k] * Z[i] + g_prev[k] * Q[i] + α_prev[k] * R[i])

    # Term 4: cross-coag mass loss, Σ_{i=1}^{k-1} (S2+yS1+y2S0) using k params
    @einsum term4_gM[k] := mask_lt1[k, i] * (d_coeff[k] * M[i] + e_coeff[k] * Z[i] + g_coeff[k] * Q[i] + α[k] * R[i])

    # Lower sums for Term 5
    @einsum M_lower[k] := mask_lt1[k, i] * M[i]
    @einsum Z_lower[k] := mask_lt1[k, i] * Z[i]

    # Assemble dM
    @einsum dM[k] = C * (N_prev[k] * Z_prev[k] + M_prev[k]^2) + C * term2_gM[k] - C * (N[k] * Z[k] + M[k]^2) - C * term4_gM[k] + C * (M[k] * M_lower[k] + N[k] * Z_lower[k]) - C * (Z[k] * N_upper[k] + M[k] * M_upper[k])

    return nothing
end
