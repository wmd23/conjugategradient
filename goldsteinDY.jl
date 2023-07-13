#
# function goldDY(x, f, gf; ϵ = 1.e-5, ϵ1 = 0.01, maxiter=400)
#----------------------------------------------------------------

# function goldDY implements the DY method (conjugate gradient) where the steplength holds the Goldstein conditions.

#
# Reference: Dai, Y.H. and Yuan, Y.X., 1998, Some properties of a new conjugate gradient method. In:Y.Yuan (Ed.)
# Advances in Nonlinear Programming (Boston: Kluwer), pp. 252-262.
#

# Input Parameters
# =================
# x (vector) vector containing the current estimation to be minimizer.
# f (function) objective function.
# gf (function) maps a point x into its gradient.

# Opcional Input Parameters
# =========================
# ϵ (Float64) tolerance given (default = 1.e-5).
# ϵ1 (Float64) constant used in the evaluation process that will decide the direction taken. (default = 0.01).
# maxiter (Int64) Set the maximum number of iterates (default = 400).

# Variables
# ==========
# ierror (Int64) Inform the following mensages: 0 - Ok! 1 - the number of iterates exceed the maximum amount allowed.
# iter (Float64) number of iterates.
# gn (Int64) number of times that the conjugate direction was not used.
# fn (Int64) number of function evaluations.
# gradfxp (vector) the gradient vector on the previous iteration.
# d (vector) descent direction.
# X (vector) contains the first coordinates of every point of the iterative process.
# Y (vector) contains the second coordinates of every point of the iterative process.
# Z (vector) contains the values of the functions avaliated in every point of the iteration process.
# gradfx (vector) the gradient vector of the current iteration.
# normgfx (Float64) norm of the gradfx.
# fx (Float64) contain the value of objective function avaliated in the current estimation.
# β (Float64) value used to adjust the descent direction.
# α (Float64) contain the value of the steplength candidate.
# iter_int (Int) number of iterations used to computate the steplength.
# x (vector) vector containing the current estimation to be minimizer.

# Output
# ======
# x (vector) vector containing the current estimation to be minimizer.
# fx (Float64) contain the value of objective function avaliated in the current estimation.
# normgfx (Float64) contain the value of the euclidean norm of the current estimation.
# iter (Int) number of iterations.
# ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
# gn (Int) number of the times where the gradient method was chosen.
# fn (Int) number of function evaluations.
# X (vector) contains the first coordinates in each estimation.
# Y (vector) contains the second coordinates in each estimation.
# Z (vector) contains the value of the function evaluation in each estimation.

include("goldstein.jl")
using LinearAlgebra

function goldDY(x, f, gf; ϵ = 1.e-5, ϵ1 = 0.01, maxiter=400)

    # initializing the variables
    ierror, iter, gn, fn = 0, 0, 0, 0            
    gradfxp, d = 0.0, 0.0    
    X, Y, Z = Float64[], Float64[], Float64[]
    push!(X, x[1]); push!(Y, x[2]); push!(Z, f(x))                   

    while true
    
        gradfx = gf(x)                          

        normgfx = norm(gradfx,2)

        if normgfx < ϵ                                              # stop condition
            fx = f(x)
            return(x,fx,normgfx,iter,ierror,gn,fn,X,Y,Z)
        end

        if iter > 0
            y = gradfx - gradfxp                                    # auxiliar computation                       
            dirtest = (d' * y)[1] - ϵ1 * norm(d,2) * normgfx        # auxiliar computation
            if dirtest < 0.0
                d = - gradfx                                        # updating the descent direction
                gn += 1
            else
                β = norm(gradfx)^2 / (d' * y)[1]                       
                d = -gradfx + β * d
            end 
        else
            d = -gradfx                                             # updating the descent direction
            gn += 1
        end

        α, ierror, iter_int = goldstein(f, gf, x, d)
        fn += iter_int                          

        x = x + α * d                                               # updating the estimation
        push!(X, x[1]); push!(Y, x[2]); push!(Z, f(x))

        iter += 1

        gradfxp = gradfx
        

        if iter > maxiter                                           # stop condition
            ierror = 1
            fx = f(x)
            return(x,fx,normgfx,iter,ierror,gn,fn,X,Y,Z)
        end

    end
end