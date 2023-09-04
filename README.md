# Main.jl
This file contains the main file which includes all the other files presented in this readme. Its code is written below.

```julia
include("armijoPRP.jl"); include("cautiousDY.jl"); include("goldsteinPRP.jl"); include("goldsteinDY.jl"); include("testfunction.jl")
 
f = powell                  # see testfunction.jl for more details 
gf = gradpowell

ndim = 100                  # set the dimension

x0 = rand(ndim)             # set the guess

x1, fx1, normx1, iter1, t1, ierror1, serror1, gn1, fn1, X1, Y1, Z1 = cautious(x0, f, gf, maxiter = 500000, ϵ1 = 0.01, ϵ = 1.0e-5)
x2, fx2, normx2, iter2, t2, ierror2, serror2, gn2, fn2, X2, Y2, Z2 = goldDY(x0, f, gf, maxiter = 500000, ϵ1 = 0.01, ϵ = 1.0e-5)
x3, fx3, normx3, iter3, t3, ierror3, serror3, gn3, fn3, X3, Y3, Z3 = armijoPRP(x0, f, gf, n, maxk = 500000, method = 0, ϵ = 1.0e-5)
x4, fx4, normx4, iter4, t4, ierror4, serror4, gn4, fn4, X4, Y4, Z4 = goldsteinPRP(x0, f, gf, n, maxk = 500000, method = 0, ϵ = 1.0e-5)
```
You only need to include these five files (armijoPRP.jl, cautiousDY.jl, goldsteinPRP.jl, goldsteinDY.jl and testfunction.jl), the first four files are algorithms to minimize $f:\mathbb{R}^n\to\mathbb{R}$ and the last one has the objective function. The variable f it is the function and gf its gradient, with ndim you can determine the dimension, $x0$ has the guess point. Moreover, the last four lines of the code show how to invoke each algorithm.

# armijoPRP.jl
This file contains the function called armijo which implements the PRP method (conjugate gradient) or the gradient method where the steplength is computated by Armijo's rule. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- ∇f (Float64) the value of the gradient avaliated in the current estimation.
- n (Int) dimension of the objective function.

Moreover, you can specify the maximum number of iterates allowed, the tolerance given, and the value of the variable called method, which is 1 by default and with this the function will implement the PRP method, otherwise the steepest descent method.

This function returns the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- fx (Float64) contain the value of objective function avaliated in the current estimation.
- normx (Float64) contain the value of the euclidean norm of the current estimation.
- iter (Int) number of iterations.
- time (Float64) the time elapsed in order to solve the minimization problem.
- ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
- serror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - at least on one iterate the steplenght is lesser than the minimum expected (steplength too small).
- counter (Int) number of the times where the gradient method was chosen.
- fn (Int) number of function evaluations.
- X (vector) contains the first coordinates in each estimation.
- Y (vector) contains the second coordinates in each estimation.
- Z (vector) contains the value of the function evaluation in each estimation.

## Remarks:

 - If ierror = 0, then a minimum was found.
 - If ndim = 2, you can use X, Y, Z to plot the sequences on level curves, for example.

# goldsteinPRP.jl
This file contains the same funcionalities of armijoPRP.jl but uses the Goldstein's line search instead.

# golsdsteinDY.jl
This file contains the function called goldDY that implements the DY method (conjugate gradient) where the steplength holds the Goldstein conditions. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- gf (function) the gradient of the objective function.

Moreover, you can specify the tolerance the constant used on the evaluation process that will decide the direction taken and the maximum number of iterations allowed.

This function returns the same information of the file armijoPRP.jl

# cautiousDY.jl
This file contains the same funcionalities of goldsteinDY.jl but uses the Armijo's line search instead.

# Armijo.jl
This file contains the function called armijo that calcutes a steplenght that Armijo's rule holds. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- gradfx (Float64) the value of the gradient avaliated in the current estimation.
- d (vector) vector containing a descent direction from the current estimation.

Moreover, you can specify the inicial steplenght $t$, the contant used on the backtracking process and the constant $c_1$ used on the evaluation of Armijo's rule $f(x+t\cdot d)\le f(x) + t\cdot c1\cdot \nabla f(x)^T\cdot d$.

This function returns the following informations:

- t (Float64) steplength which satisfies the Armijo's rule.
- iter (Int) number of iterations.

## Remarks:

- As you might know, you can use as many optional parameters to invoke the function as you want.
- $t > 0$.
- $y, c1 \in (0, 1)$.

# Goldstein.jl
This file contains the function called goldstein that implements the goldstein line search. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- ∇f (function) the gradient of the objective function.
- d (vector) vector containing a descent direction from the current estimation.

Moreover, you can specify the initial steplenght the constant used on the evaluation of Goldstein conditions and maximum number of iterations allowed.

This function returns the following informations:

- α (Float64) steplength which satisfies the goldstein conditions.
- ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
- fn (Int) number of function evaluations.
