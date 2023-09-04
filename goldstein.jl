#
# Adaptation of strongwolfe.jl
#

# function goldstein(f, ∇, x, d; α=1, β=1e-4, σ=0.1, maxiter = 1000)
#-----------------------------------------------
# function goldstein implements the goldstein line search


# Input Parameters
# =================
# f (function) objective function.
# ∇ (function) the gradient of the objective function.
# x (vector) vector containing the current estimation to be minimizer.
# d (vector) vector containing a descent direction from the current estimation.

# Optional Input Parameters
# =========================
# α (Float64) inicial steplength. (default α=1)
# β (Float64) used on the evaluation of Goldstein conditions. (default β=1e-4).
# σ (Float64) used on the evaluation of Goldstein conditions. (default σ=0.1).
# maxiter (Int) maximum number of iterations allowed. (default = 1000).

# Variables
# ==========
# y0 (Float64) y0 = φ(0), where φ(α) = f(x + α*d).
# g0 (Float64) g0 = φ'(0).
# y_prev (Float64) the value of φ on the previous iteration.
# α_prev (Float64) the previous estimative for steplength.
# αlo (Float64) lesser endpoint of an interval which contain steplengths where Goldstein conditions holds.
# αhi (Float64) greater endpoint of an interval which contain steplengths where Goldstein conditions holds.
# ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
# iter (Int) number of iterations used to computate the steplength.
# y (Float64) y = φ(α)
# fn (Int) number of function evaluations.
# g (Float64) g = φ'(α)
# α (Float64) contain the value of the steplength candidate.
# ylo (Float64) ylo = φ(αlo)

# Output
# ======
# α (Float64) steplength which satisfies the goldstein conditions.
# ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
# fn (Int) number of function evaluations.

using LinearAlgebra

function goldstein(f, ∇, x, d; α=1, β=1e-4, σ=0.1, maxiter = 1000)
    y0, g0, y_prev, α_prev = f(x), ∇(x)⋅d, NaN, 0
    αlo, αhi = NaN, NaN
    ierror, iter, fn = 0, 0, 1

    #bracket phase
    while true
        y = f(x + α*d)
        fn += 1
        if y0 + (1-β)*α*g0 > y || (!isnan(y_prev) && y ≥ y_prev)
            αlo, αhi = α_prev, α
            break
        end
        g = ∇(x + α*d)⋅d
        if y ≤ y0 + β*α*g0
            return α, ierror, fn
        elseif g ≥ 0
            αlo, αhi = α, α_prev
            break
        end
        y_prev, α_prev, α = y, α, 2α
    end

    # zoom phase
    ylo = f(x + αlo*d)
    fn += 1
    while true
        α = (αlo + αhi)/2
        iter = iter + 1
        if iter > maxiter
            ierror = 1
            return α, ierror, fn
        end
        y = f(x + α*d)
        fn += 1
        if y0 + (1-β)*α*g0 > y || y0 + (1-β)*α*g0 ≥ ylo
            αhi = α
        else
            g = ∇(x + α*d)⋅d
            if y ≤ y0 + β*α*g0
                return α, ierror, fn

            elseif g*(αhi - αlo) ≥ 0
                αhi = αlo
            end
            αlo = α
        end
    end
end