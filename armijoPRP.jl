# function armijoPRP(x, f, ∇f, n; maxk = 1000, ϵ = 1e-5, method = 1)
#-----------------------------------------------
# function armijoPRP implements the PRP method (conjugate gradient) or the gradient method where the steplength is computated by Armijo's rule

# Input Parameters
# =================
# x (vector) vector containing the current estimation to be minimizer.
# f (function) objective function.
# ∇f (function) maps a point x into its gradient.
# n (Int) dimension of the objective function.

# Optional Input Parameters
# =========================
# maxk (Int) set the maximum number of iterations. (default = 1000).
# ϵ (Float64) tolerance given. (default = 1e-5).
# method (Int) if method = 1 the function will implement the PRP method, otherwise the gradient method.

# Variables
# ==========
# d (vector) vector containing a descent direction from the current estimation.
# gradfx (Float64) the value of the gradient avaliated in the previous estimation.
# ierror (Int) the value stored in this variable informs the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
# k (Int) number of iterations.
# fn (Int) number of function evaluations.
# counter (Int) number of the times where the gradient method was chosen.
# X (vector) contains the first coordinates in each estimation.
# Y (vector) contains the second coordinates in each estimation.
# Z (vector) contains the value of the function evaluation in each estimation.
# t (Float64) contain the value of the steplength candidate.
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

include("armijo.jl")

using LinearAlgebra

function armijoPRP(x, f, ∇f, n; maxk = 1000, ϵ = 1e-5, method = 1)

    # initializing the variables
    d = -∇f(x)                  
    gradfx = -d
    ierror = 0
    k = 0
    fn = 0
    counter = 1;
    X, Y, Z = Float64[], Float64[], Float64[]

    push!(X, x[1]); push!(Y, x[2]); push!(Z, f(x))  # adding the correponding values to each vector

    while norm(∇f(x)) > ϵ && k < maxk               # stop conditions  

        t, iter = armijo(x, f, gradfx, d)
        fn += iter + 1                              # the exact number of the functions evaluations

        x = x + t*d                                 # updating the estimation
        push!(X, x[1]); push!(Y, x[2]); push!(Z, f(x))

        newgradfx = ∇f(x)
        
        
        if method == 1
            # restart
            if (k+1)%n != 0 
                βk = (newgradfx'⋅(newgradfx - gradfx))/(gradfx'⋅gradfx)
                βk = max(0, βk)                                         # because βk can be negative
                if βk == 0
                    counter = counter + 1
                end
            else
                βk = 0
                counter = counter + 1
            end
        else
            # gradient method
            βk = 0
            counter = counter + 1
        end

        d = -newgradfx + βk*d           # updating the descent direction
        k = k + 1
        gradfx = newgradfx
    end

    if k < maxk
        # gathering up the last two useful informations
        fx = f(x)
        normx = norm(gradfx)
        return x, fx, normx, k, ierror, counter, fn, X, Y, Z
    else
        # gathering up the last three useful informations
        ierror = 1
        fx = f(x)
        normx = norm(gradfx)
        return x, fx, normx, k, ierror, counter, fn, X, Y, Z
    end
end