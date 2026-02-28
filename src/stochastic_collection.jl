export StochasticCollectionCoalescence

using SymbolicUtils: @arrayop, @makearray

# Helper functions for @arrayop — take scalar symbolic arguments, return symbolic expressions.

function _sce_Z(n_dl, m_dl, xi_p)
    return ifelse(n_dl > 0, xi_p * m_dl^2 / n_dl, 0.0)
end

function _sce_Q(n_dl, m_dl, xi_p)
    return ifelse(n_dl > 0, xi_p^2 * m_dl^3 / n_dl^2, 0.0)
end

function _sce_R(n_dl, m_dl, xi_p)
    return ifelse(n_dl > 0, xi_p^3 * m_dl^4 / n_dl^3, 0.0)
end

function _sce_f(n_dl, m_dl, xk)  # Eq. 13a + positivity Eq. 15
    xb = ifelse(n_dl > 0, m_dl / n_dl, xk)
    return ifelse(
        xb < xk, 0.0,
        ifelse(
            xb > 2 * xk, 2 * n_dl / xk,
            4 * n_dl / xk - 2 * m_dl / xk^2
        )
    )
end

function _sce_psi(n_dl, m_dl, xk)  # Eq. 13b + positivity Eq. 15
    xb = ifelse(n_dl > 0, m_dl / n_dl, xk)
    return ifelse(
        xb < xk, 2 * n_dl / xk,
        ifelse(
            xb > 2 * xk, 0.0,
            2 * m_dl / xk^2 - 2 * n_dl / xk
        )
    )
end

function _sce_alpha(n_dl, m_dl, xk)
    return (_sce_f(n_dl, m_dl, xk) - _sce_psi(n_dl, m_dl, xk)) / (2 * xk)
end

"""
    StochasticCollectionCoalescence(; name=:StochasticCollectionCoalescence, I=36, x1=1.6e-14, kernel_type=:constant, kernel_params=Dict(:K0 => 1e-10))

Solve the stochastic collection equation (SCE) for cloud droplet coalescence using
a two-moment method in discrete mass categories with `p = 2` (mass-doubling categories).

**Reference**: Tzivion, S., Feingold, G., and Levin, Z. (1987) "A Numerical Solution of the Kinetic Collection Equation", *Journal of the Atmospheric Sciences*, 44(24), 3827-3844; and Tzivion, S., Feingold, G., and Levin, Z. (1989) "The Evolution of Raindrop Spectra. Part II: Collisional Collection/Breakup and Evaporation in a Rainshaft", *Journal of the Atmospheric Sciences*, 46(21), 3312-3327.

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
SymbolicUtils.jl) for all computations: closure relations, linear distribution parameters,
shifted arrays, and summation reductions. Rate assembly uses level-2 `@arrayop` with
`ifelse` masks for variable-bound sums and the `ifelse(i==1, ...)` trick for per-k terms.

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
        name = :StochasticCollectionCoalescence,
        I::Int = 36, x1::Float64 = 1.6e-14,
        kernel_type::Symbol = :constant,
        kernel_params::Dict = Dict(:K0 => 1.0e-10)
    )

    # ---- Validate kernel ----
    if kernel_type == :constant
        K0_val = Float64(kernel_params[:K0])
    elseif kernel_type == :golovin
        C_val = Float64(kernel_params[:C])
    else
        error("Unsupported kernel_type: $kernel_type. Use :constant or :golovin.")
    end

    @constants begin
        XI_P = 1.0625, [description = "Closure parameter ξ_p for relating higher-order moments (dimensionless)", unit = u"1"]
        p = 2.0, [description = "Mass doubling factor for geometric category spacing (dimensionless)", unit = u"1"]
        one_m3 = 1.0, [description = "Unit volume for dimensional analysis", unit = u"m^3"]
        one_kg = 1.0, [description = "Unit mass for dimensional analysis", unit = u"kg"]
        one_s = 1.0, [description = "Unit time for dimensional analysis", unit = u"s"]
    end

    @variables begin
        Nk(t)[1:I], [description = "Number concentration in category k", unit = u"m^-3"]
        Mk(t)[1:I], [description = "Mass concentration in category k", unit = u"kg*m^-3"]
    end

    # ---- Closure relations via @arrayop (level 1) ----
    Z_arr = @arrayop (k,) _sce_Z(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), XI_P) k in 1:I
    Q_arr = @arrayop (k,) _sce_Q(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), XI_P) k in 1:I

    # ---- Linear distribution parameters via @arrayop (level 1, Eq. 13a,b + Eq. 15a,b) ----
    f_arr = @arrayop (k,) _sce_f(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), x1 * p^(k - 1)) k in 1:I
    ψ_arr = @arrayop (k,) _sce_psi(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), x1 * p^(k - 1)) k in 1:I
    α_arr = @arrayop (k,) _sce_alpha(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), x1 * p^(k - 1)) k in 1:I

    # ---- Shifted arrays via @makearray (level 1) ----
    N_prev_arr = @makearray N_prev_arr[1:I] begin
        N_prev_arr[1:1] => [0.0]
        N_prev_arr[2:I] => @arrayop (k,) Nk[k] * one_m3 k in 1:(I - 1)
    end
    M_prev_arr = @makearray M_prev_arr[1:I] begin
        M_prev_arr[1:1] => [0.0]
        M_prev_arr[2:I] => @arrayop (k,) Mk[k] * (one_m3 / one_kg) k in 1:(I - 1)
    end

    ψ_prev_arr = @makearray ψ_prev_arr[1:I] begin
        ψ_prev_arr[1:1] => [0.0]
        ψ_prev_arr[2:I] => @arrayop (k,) _sce_psi(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), x1 * p^(k - 1)) k in 1:(I - 1)
    end
    α_prev_arr = @makearray α_prev_arr[1:I] begin
        α_prev_arr[1:1] => [0.0]
        α_prev_arr[2:I] => @arrayop (k,) _sce_alpha(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), x1 * p^(k - 1)) k in 1:(I - 1)
    end
    f_prev_arr = @makearray f_prev_arr[1:I] begin
        f_prev_arr[1:1] => [0.0]
        f_prev_arr[2:I] => @arrayop (k,) _sce_f(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), x1 * p^(k - 1)) k in 1:(I - 1)
    end

    # ---- Kernel-specific rate computation via level-2 @arrayop ----
    if kernel_type == :constant
        # dN: Eq. 9a with K = K0
        dN_arr = @arrayop (k,) (
            # Per-k terms (Terms 1+3), gated to fire once via i==1
            ifelse(
                i == 1,
                K0_val / 2 * N_prev_arr[k]^2 - K0_val * (Nk[k] * one_m3)^2,
                0.0
            ) +
                # Term 2: sum_{i=1}^{k-2} with k-1 shifted params
                ifelse(
                i <= k - 2,
                K0_val * (ψ_prev_arr[k] * (Mk[i] * (one_m3 / one_kg)) + α_prev_arr[k] * Z_arr[i]),
                0.0
            ) +
                # Term 4: -sum_{i=1}^{k-1}
                ifelse(
                i <= k - 1,
                -K0_val * (ψ_arr[k] * (Mk[i] * (one_m3 / one_kg)) + α_arr[k] * Z_arr[i]),
                0.0
            ) +
                # Folded -N_dl[k]*N_upper[k]
                ifelse(
                i >= k + 1,
                -K0_val * (Nk[k] * one_m3) * (Nk[i] * one_m3),
                0.0
            )
        ) k in 1:I i in 1:I

        # dM: Eq. 9b with K = K0
        dM_arr = @arrayop (k,) (
            # Per-k terms (Terms 1+3)
            ifelse(
                i == 1,
                K0_val * N_prev_arr[k] * M_prev_arr[k]
                    - K0_val * (Nk[k] * one_m3) * (Mk[k] * (one_m3 / one_kg)),
                0.0
            ) +
                # Term 2: sum_{i=1}^{k-2} with shifted coefficients inlined
                ifelse(
                i <= k - 2,
                K0_val * (
                    2 * x1 * p^(k - 2) * ψ_prev_arr[k] * (Mk[i] * (one_m3 / one_kg))
                        + f_prev_arr[k] / 2 * Z_arr[i]
                        + α_prev_arr[k] * Q_arr[i]
                ),
                0.0
            ) +
                # Terms 4+5 combined (same mask i<=k-1): -sum + folded M_lower
                ifelse(
                i <= k - 1,
                -K0_val * (
                    2 * x1 * p^(k - 1) * ψ_arr[k] * (Mk[i] * (one_m3 / one_kg))
                        + f_arr[k] / 2 * Z_arr[i]
                        + α_arr[k] * Q_arr[i]
                )
                    + K0_val * (Nk[k] * one_m3) * (Mk[i] * (one_m3 / one_kg)),
                0.0
            ) +
                # Folded -M_dl[k]*N_upper[k]
                ifelse(
                i >= k + 1,
                -K0_val * (Mk[k] * (one_m3 / one_kg)) * (Nk[i] * one_m3),
                0.0
            )
        ) k in 1:I i in 1:I

    else  # :golovin
        R_arr = @arrayop (k,) _sce_R(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), XI_P) k in 1:I

        Z_prev_arr = @makearray Z_prev_arr[1:I] begin
            Z_prev_arr[1:1] => [0.0]
            Z_prev_arr[2:I] => @arrayop (k,) _sce_Z(Nk[k] * one_m3, Mk[k] * (one_m3 / one_kg), XI_P) k in 1:(I - 1)
        end

        # dN: Eq. 9a with K = C(x+y)
        dN_arr = @arrayop (k,) (
            # Per-k terms
            ifelse(
                i == 1,
                C_val * N_prev_arr[k] * M_prev_arr[k]
                    - 2 * C_val * (Nk[k] * one_m3) * (Mk[k] * (one_m3 / one_kg)),
                0.0
            ) +
                # Term 2: sum_{i=1}^{k-2} with shifted coefficients
                ifelse(
                i <= k - 2,
                C_val * (
                    2 * x1 * p^(k - 2) * ψ_prev_arr[k] * (Mk[i] * (one_m3 / one_kg))
                        + f_prev_arr[k] / 2 * Z_arr[i]
                        + α_prev_arr[k] * Q_arr[i]
                ),
                0.0
            ) +
                # Term 4: -sum_{i=1}^{k-1}
                ifelse(
                i <= k - 1,
                -C_val * (
                    2 * x1 * p^(k - 1) * ψ_arr[k] * (Mk[i] * (one_m3 / one_kg))
                        + f_arr[k] / 2 * Z_arr[i]
                        + α_arr[k] * Q_arr[i]
                ),
                0.0
            ) +
                # Folded -M_dl[k]*N_upper[k] - N_dl[k]*M_upper[k]
                ifelse(
                i >= k + 1,
                -C_val * (
                    (Mk[k] * (one_m3 / one_kg)) * (Nk[i] * one_m3)
                        + (Nk[k] * one_m3) * (Mk[i] * (one_m3 / one_kg))
                ),
                0.0
            )
        ) k in 1:I i in 1:I

        # dM: Eq. 9b with K = C(x+y)
        dM_arr = @arrayop (k,) (
            # Per-k terms
            ifelse(
                i == 1,
                C_val * (N_prev_arr[k] * Z_prev_arr[k] + M_prev_arr[k]^2)
                    - C_val * ((Nk[k] * one_m3) * Z_arr[k] + (Mk[k] * (one_m3 / one_kg))^2),
                0.0
            ) +
                # Term 2: sum_{i=1}^{k-2} with shifted coefficients
                ifelse(
                i <= k - 2,
                C_val * (
                    4 * (x1 * p^(k - 2))^2 * ψ_prev_arr[k] * (Mk[i] * (one_m3 / one_kg))
                        + x1 * p^(k - 2) * (f_prev_arr[k] / 2 + 2 * ψ_prev_arr[k]) * Z_arr[i]
                        + (f_prev_arr[k] - ψ_prev_arr[k]) * Q_arr[i]
                        + α_prev_arr[k] * R_arr[i]
                ),
                0.0
            ) +
                # Terms 4+5+6 combined (mask i<=k-1): -term4 + M_lower + Z_lower folded
                ifelse(
                i <= k - 1,
                -C_val * (
                    4 * (x1 * p^(k - 1))^2 * ψ_arr[k] * (Mk[i] * (one_m3 / one_kg))
                        + x1 * p^(k - 1) * (f_arr[k] / 2 + 2 * ψ_arr[k]) * Z_arr[i]
                        + (f_arr[k] - ψ_arr[k]) * Q_arr[i]
                        + α_arr[k] * R_arr[i]
                )
                    + C_val * (Mk[k] * (one_m3 / one_kg)) * (Mk[i] * (one_m3 / one_kg))
                    + C_val * (Nk[k] * one_m3) * Z_arr[i],
                0.0
            ) +
                # Folded -Z[k]*N_upper[k] - M_dl[k]*M_upper[k]
                ifelse(
                i >= k + 1,
                -C_val * (
                    Z_arr[k] * (Nk[i] * one_m3)
                        + (Mk[k] * (one_m3 / one_kg)) * (Mk[i] * (one_m3 / one_kg))
                ),
                0.0
            )
        ) k in 1:I i in 1:I
    end

    # ---- Build equations with proper units (using symbolic array indexing, no scalarize) ----
    eqs = vcat(
        [D(Nk[k]) ~ dN_arr[k] / (one_m3 * one_s) for k in 1:I],               # m⁻³ s⁻¹
        [D(Mk[k]) ~ dM_arr[k] * one_kg / (one_m3 * one_s) for k in 1:I],  # kg m⁻³ s⁻¹
    )

    return System(eqs, t; name, checks = ModelingToolkit.CheckComponents)
end
