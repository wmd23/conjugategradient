# function armijo(x, f, gradfx, d; t = 0.5, y=0.5, c1=0.1)
#---------------------------------------------------------
# function armijo implements the backtracking algorithm. (Armijo)

#
# Reference: NOCEDAL, Jorge; WRIGHT, Stephen J. (Ed.). Numerical optimization. New York, NY: Springer New York, 1999.
#

# Input Parameters
# =================
# x (vector) containing the current estimation to be minimizer.
# f (function) objective function.
# gradfx (Float64) the value of the gradient avaliated in the current estimation.
# d (vector) vector containing a descent direction from the current estimation.

# Optional Input Parameters
# =========================
# t (Float64) inicial steplength. (default = 0.5).
# y (Float64) constant used on the backtracking process. (default = 0.5).
# c1 (Float64) Used on the evaluation of Armijo's rule (f(x+t*d) ≤ fx + t*c1*gradfx'*d) (default = 0.1).

# Variables
# ==========
# iter (Int) number of iterations.
# fx (Float64) contain the value of objective function avaliated in the current estimation.
# gd (Float64) contain the value of the value of φ'(0) where φ(t) = f(x + t*d).
# t (Float64) contain the value of the steplength candidate.

# Output
# ======
# t (Float64) steplength which satisfies the Armijo's rule.
# iter (Int) number of iterations.

# Remark
# =======
# Use this function only when the global convergence of the method is assured, otherwise an infinite loop
# can be achieved
#

using LinearAlgebra

function armijo(x, f, gradfx, d; t = 0.5, y=0.5, c1=0.1)
    
    iter = 0                        # set the counter
    fx = f(x)                       # this prevents unecessary function evaluations
    gd = gradfx'*d; gd = gd[1]      # auxiliar computations

    while f(x+t*d) > fx + t*c1*gd   # while the Armijo's rule is not hold
        t = y*t                     # backtracking process
        iter += 1                   # updating the number of iterations
    end

    return t, iter                  # return the steplength that holds Armijo's rule

end