export StochasticCollectionCoalescence

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
    eqs = Equation[]
    for k in 1:I
        push!(eqs, D(Nk[k]) ~ result[k] / (one_m3 * one_s))                    # m⁻³ s⁻¹
        push!(eqs, D(Mk[k]) ~ result[I + k] * one_kg / (one_m3 * one_s))       # kg m⁻³ s⁻¹
    end

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
# Implements Eq. (9a,b) with K constant.
# ============================================================================
function _sce_constant_rhs!(dN, dM, N, M, I, xb, K0)
    @inbounds for k in 1:I
        xk = xb[k]
        dNk = 0.0
        dMk = 0.0

        # ---- dN_k/dt (Eq. 9a) ----

        # Term 1: Autoconversion gain from k-1 (coag of two k-1 particles → k)
        if k >= 2
            dNk += (K0 / 2) * N[k-1]^2
        end

        # Term 2: Cross-coagulation gain (particles from i < k-1 combining with k-1)
        if k >= 3
            xkm1 = xb[k-1]
            fkm1, ψkm1 = _sce_linear_params(N[k-1], M[k-1], xkm1)
            for i in 1:(k-2)
                Zi = _sce_Z(M[i], N[i])
                dNk += K0 * _sce_S0(fkm1, ψkm1, xkm1, M[i], Zi)
            end
        end

        # Term 3: Autoconversion loss from k (coag of two k particles → k+1)
        dNk -= (K0 / 2) * N[k]^2

        # Term 4: Cross-coagulation loss (particles from i < k combining with k)
        if k >= 2
            fk, ψk = _sce_linear_params(N[k], M[k], xk)
            for i in 1:(k-1)
                Zi = _sce_Z(M[i], N[i])
                dNk -= K0 * _sce_S0(fk, ψk, xk, M[i], Zi)
            end
        end

        # Term 5: Self-coagulation number loss (k with k, full integral)
        # Factor 1/2 because ∫∫ n_k(y) n_k(x) K dx dy double-counts pairs
        dNk -= (K0 / 2) * N[k]^2

        # Term 6: Number loss from coagulation with higher categories
        for i in (k+1):I
            dNk -= K0 * N[k] * N[i]
        end

        dN[k] = dNk

        # ---- dM_k/dt (Eq. 9b) ----

        # Term 1: Autoconversion mass gain from k-1
        if k >= 2
            dMk += K0 * N[k-1] * M[k-1]
        end

        # Term 2: Cross-coagulation mass gain (incomplete integrals over k-1)
        if k >= 3
            xkm1 = xb[k-1]
            fkm1, ψkm1 = _sce_linear_params(N[k-1], M[k-1], xkm1)
            for i in 1:(k-2)
                Zi = _sce_Z(M[i], N[i])
                Qi = _sce_Q(M[i], N[i])
                dMk += K0 * (_sce_S1(fkm1, ψkm1, xkm1, M[i], Zi) +
                             _sce_yS0(fkm1, ψkm1, xkm1, Zi, Qi))
            end
        end

        # Term 3: Autoconversion mass loss from k
        dMk -= K0 * N[k] * M[k]

        # Term 4: Cross-coagulation mass loss (incomplete integrals over k)
        if k >= 2
            fk, ψk = _sce_linear_params(N[k], M[k], xk)
            for i in 1:(k-1)
                Zi = _sce_Z(M[i], N[i])
                Qi = _sce_Q(M[i], N[i])
                dMk -= K0 * (_sce_S1(fk, ψk, xk, M[i], Zi) +
                             _sce_yS0(fk, ψk, xk, Zi, Qi))
            end
        end

        # Term 5: Mass gain from coagulation with lower categories
        for i in 1:(k-1)
            dMk += K0 * N[k] * M[i]
        end

        # Term 6: Mass loss from coagulation with higher categories
        for i in (k+1):I
            dMk -= K0 * M[k] * N[i]
        end

        dM[k] = dMk
    end
    return nothing
end

# ============================================================================
# Numerical RHS for Golovin kernel K(x,y) = C*(x+y)
# Implements Eq. (9a,b) with K = C(x+y).
# ============================================================================
function _sce_golovin_rhs!(dN, dM, N, M, I, xb, C)
    @inbounds for k in 1:I
        xk = xb[k]
        dNk = 0.0
        dMk = 0.0

        # ---- dN_k/dt (Eq. 9a with K=C(x+y)) ----

        # Term 1: autoconversion gain from k-1
        if k >= 2
            dNk += C * N[k-1] * M[k-1]
        end

        # Term 2: cross-coag gain (incomplete integrals over k-1)
        if k >= 3
            xkm1 = xb[k-1]
            fkm1, ψkm1 = _sce_linear_params(N[k-1], M[k-1], xkm1)
            for i in 1:(k-2)
                Zi = _sce_Z(M[i], N[i])
                Qi = _sce_Q(M[i], N[i])
                dNk += C * (_sce_S1(fkm1, ψkm1, xkm1, M[i], Zi) +
                            _sce_yS0(fkm1, ψkm1, xkm1, Zi, Qi))
            end
        end

        # Term 3: autoconversion loss
        dNk -= C * N[k] * M[k]

        # Term 4: cross-coag loss (incomplete integrals over k)
        if k >= 2
            fk, ψk = _sce_linear_params(N[k], M[k], xk)
            for i in 1:(k-1)
                Zi = _sce_Z(M[i], N[i])
                Qi = _sce_Q(M[i], N[i])
                dNk -= C * (_sce_S1(fk, ψk, xk, M[i], Zi) +
                            _sce_yS0(fk, ψk, xk, Zi, Qi))
            end
        end

        # Term 5: self-coag number loss (factor 1/2 for pair double-counting)
        # (1/2) C ∫∫ n_k(x) n_k(y) (x+y) dx dy = (1/2) C * 2 N_k M_k = C N_k M_k
        dNk -= C * N[k] * M[k]

        # Term 6: number loss from coag with higher categories
        for i in (k+1):I
            dNk -= C * (N[i] * M[k] + M[i] * N[k])
        end

        dN[k] = dNk

        # ---- dM_k/dt (Eq. 9b with K=C(x+y)) ----

        Zk = _sce_Z(M[k], N[k])

        # Term 1: autoconversion mass gain
        if k >= 2
            Zkm1 = _sce_Z(M[k-1], N[k-1])
            dMk += C * (N[k-1] * Zkm1 + M[k-1]^2)
        end

        # Term 2: cross-coag mass gain (incomplete integrals over k-1)
        if k >= 3
            xkm1 = xb[k-1]
            fkm1, ψkm1 = _sce_linear_params(N[k-1], M[k-1], xkm1)
            for i in 1:(k-2)
                Zi = _sce_Z(M[i], N[i])
                Qi = _sce_Q(M[i], N[i])
                Ri = _sce_R(M[i], N[i])
                S2 = 4 * xkm1^2 * ψkm1 * M[i] + (xkm1 * fkm1 - 4 * xkm1 * ψkm1) / 2 * Zi
                yS1 = 2 * (2 * xkm1 * ψkm1 * Zi + (fkm1 - 2 * ψkm1) / 2 * Qi)
                y2S0 = ψkm1 * Qi + (fkm1 - ψkm1) / (2 * xkm1) * Ri
                dMk += C * (S2 + yS1 + y2S0)
            end
        end

        # Term 3: autoconversion mass loss
        dMk -= C * (N[k] * Zk + M[k]^2)

        # Term 4: cross-coag mass loss (incomplete integrals over k)
        if k >= 2
            fk, ψk = _sce_linear_params(N[k], M[k], xk)
            for i in 1:(k-1)
                Zi = _sce_Z(M[i], N[i])
                Qi = _sce_Q(M[i], N[i])
                Ri = _sce_R(M[i], N[i])
                S2 = 4 * xk^2 * ψk * M[i] + (xk * fk - 4 * xk * ψk) / 2 * Zi
                yS1 = 2 * (2 * xk * ψk * Zi + (fk - 2 * ψk) / 2 * Qi)
                y2S0 = ψk * Qi + (fk - ψk) / (2 * xk) * Ri
                dMk -= C * (S2 + yS1 + y2S0)
            end
        end

        # Term 5: mass gain from coag with lower categories
        for i in 1:(k-1)
            dMk += C * (M[i] * M[k] + _sce_Z(M[i], N[i]) * N[k])
        end

        # Term 6: mass loss from coag with higher categories
        for i in (k+1):I
            dMk -= C * (N[i] * Zk + M[i] * M[k])
        end

        dM[k] = dMk
    end
    return nothing
end
