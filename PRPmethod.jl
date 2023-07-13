# function conjugadoPRP(x, f, ∇f, n; maxk, ϵ)
#-----------------------------------------------
# function conjugadoPRP implements the PRP method (conjugate gradient) where the steplength is computated by strong wolfe conditions.

# Input Parameters
# =================
# x (vector) vector containing the current estimation to be minimizer.
# f (function) objective function.
# ∇f (function) the gradient of the objective function.
# n (Int) dimension of the objective function.

# Optional Input Parameters
# =========================
# maxk (Int) set the maximum number of iterations. (default = 1000).
# ϵ (Float64) tolerance given. (default = 0.5).

# Variables
# ==========
# d (vector) vector containing a descent direction from the current estimation.
# gradfx (Float64) the value of the gradient avaliated in the previous estimation.
# ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
# k (Int) number of iterations.
# fn (Int) number of function evaluations.
# counter (Int) number of the times where the gradient method was chosen.
# X (vector) contains the first coordinates in each estimation.
# Y (vector) contains the second coordinates in each estimation.
# Z (vector) contains the value of the function evaluation in each estimation.
# t (Float) steplength computated by strong wolfe conditions conditions.
# iter (Int) number of iterations used to computate the steplength.
# x (vector) vector containing the current estimation to be minimizer.
# newgradfx (Float64) the value of the gradient avaliated in the current estimation. 
# βk (Float64) value used to adjust the descent direction.
# fx (Float64) contain the value of objective function avaliated in the current estimation.
# normx (Float64) contain the value of the euclidean norm of the current estimation.

# Output
# ======
# x (vector) vector containing the current estimation to be minimizer.
# fx (Float64) contain the value of objective function avaliated in the current estimation.
# normx (Float64) contain the value of the euclidean norm of the current estimation.
# k (Int) number of iterations.
# ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
# counter (Int) number of the times where the gradient method was chosen.
# fn (Int) number of function evaluations.
# X (vector) contains the first coordinates in each estimation.
# Y (vector) contains the second coordinates in each estimation.
# Z (vector) contains the value of the function evaluation in each estimation.

include("strongwolf.jl")
using LinearAlgebra

function conjugadoPRP(x, f, ∇f, n; maxk = 1000, ϵ = 1e-5)
    d = -∇f(x)
    gradfx = -d
    ierror = 0
    k = 0
    fn = 0
    counter = 1;
    X, Y, Z = Float64[], Float64[], Float64[]
    push!(X, x[1]); push!(Y, x[2]); push!(Z, f(x))

    while norm(∇f(x)) > ϵ && k < maxk

        t, error, fn_int = wolfe(f, ∇f, x, d)
        fn += fn_int

        x = x + t*d
        push!(X, x[1]); push!(Y, x[2]); push!(Z, f(x))

        newgradfx = ∇f(x)
        
        if (k+1)%n != 0 
            βk = (newgradfx'⋅(newgradfx - gradfx))/(gradfx'⋅gradfx)
            βk = max(0, βk)
            if βk == 0
                counter = counter + 1
            end
        else
            βk = 0
            counter = counter + 1
        end
        d = -newgradfx + βk*d
        k = k + 1
        gradfx = newgradfx
    end

    if k < maxk
        fx = f(x)
        normx = norm(gradfx)
        return x, fx, normx, k, ierror, counter, fn, X, Y, Z
    else
        ierror = 1
        fx = f(x)
        normx = norm(gradfx)
        return x, fx, normx, k, ierror, counter, fn, X, Y, Z
    end
end