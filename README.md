# Armijo.jl
This file contains the function called armijo that calcutes a steplenght that Armijo's rule holds. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- gradfx (Float64) the value of the gradient avaliated in the current estimation.
- d (vector) vector containing a descent direction from the current estimation.

Moreover, you can specify the inicial steplenght $t$, the contant used on the backtracking process and the constant $c_1$ used on the evaluation of Armijo's rule $f(x+t\cdot d)\le f(x) + t\cdot c1\cdot \nabla f(x)^T\cdot d$.

## Example:

```julia
  include("testfunctions.jl")   # see testfunctions.jl code for more details
  include("armijo.jl")          # see armijo.jl code for more details
  ndim = 100                    # dimension
  x = rand(ndim)
  f = pen_I                 # this function is in the testfunctions.jl file
  gf = gradpen_I            # this function is also in the testfunctions.jl file and gf means the gradient of f
  gradfx = gradpen_I(x)
  d = -gradfx                                                       # we can use this d as a descent direction
  t1, iter1 = armijo(x, f, gradfx, d)                               # without using optional parameters
  t2, iter2 = armijo(x, f, gradfx, d, t = 1.0, y = 0.9, c1 = 0.2)   # using optional parameters
```

This function returns the following informations:

- t (Float64) steplength which satisfies the Armijo's rule.
- iter (Int) number of iterations.

## Remarks:

- As you might know, you can use as many optional parameters to invoke the function as you want.
- $t > 0$.
- $y, c1 \in (0, 1)$.

# armijoPRP.jl
This file contains the function called conjugadoPRP that implements the PRP method (conjugate gradient) or the gradient method where the steplength is computated by Armijo's rule. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- ∇f (function) the gradient of the objective function.
- n (Int) dimension of the objective function.

Moreover, you can specify the maximum number of iterations allowed the tolerance and the method.

## Example: 

```julia
  include("testfunctions.jl")   # see testfunctions.jl code for more details
  include("armijoPRP.jl")       # see armijoPRP.jl code for more details
  ndim = 100                    # dimension
  x = rand(ndim)
  x1, fx1, normx1, iter1, ierror1, counter1, fn1, X1, Y1, Z1 = armijoPRP(x, pen_I, gradpen_I, ndim)
  x2, fx2, normx2, iter2, ierror2, counter2, fn2, X2, Y2, Z2 = armijoPRP(x, pen_I, gradpen_I, ndim, maxk = 400, ϵ = 1.e-6, method = 0)
```

This function returns the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- fx (Float64) contain the value of objective function avaliated in the current estimation.
- normx (Float64) contain the value of the euclidean norm of the current estimation.
- iter (Int) number of iterations.
- ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
- counter (Int) number of the times where the gradient method was chosen.
- fn (Int) number of function evaluations.
- X (vector) contains the first coordinates in each estimation.
- Y (vector) contains the second coordinates in each estimation.
- Z (vector) contains the value of the function evaluation in each estimation.

 ## Remarks:

 - If error = 0, then a minimum was found.
 - If ndim = 2, you can use X, Y, Z to plot the sequences on level curves, for example.
 - $maxk \in \mathbb{N}$.
 - $ϵ > 0$.
 - If method $\ne$ 1, then the gradient method will be choosed.

# cautiousDY.jl
This file contains the function called cautious that implements the cautious DY method. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- gf (function) the gradient of the objective function.

Moreover, you can specify the tolerance the constant used on the evaluation process that will decide the direction taken and the maximum number of iterations allowed.

## Example: 

```julia
  include("testfunctions.jl")   # see testfunctions.jl code for more details
  include("cautiousDY.jl")      # see cautiousDY.jl code for more details 
  ndim = 100                    # dimension
  x = rand(ndim)
  x1, fx1, normx1, iter1, ierror1, counter1, fn1, X1, Y1, Z1 = cautious(x, pen_I, gradpen_I)
  x2, fx2, normx2, iter2, ierror2, counter2, fn2, X2, Y2, Z2 = cautious(x, pen_I, gradpen_I, ϵ = 1.e-6, ϵ1 = 0.10, maxiter=500)
```

This function returns the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- fx (Float64) contain the value of objective function avaliated in the current estimation.
- normx (Float64) contain the value of the euclidean norm of the current estimation.
- iter (Int) number of iterations.
- ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
- counter (Int) number of the times where the gradient method was chosen.
- fn (Int) number of function evaluations.
- X (vector) contains the first coordinates in each estimation.
- Y (vector) contains the second coordinates in each estimation.
- Z (vector) contains the value of the function evaluation in each estimation.

## Remarks:

 - If error = 0, then a minimum was found.
 - If ndim = 2, you can use X, Y, Z to plot the sequences on level curves, for example.

# Goldstein.jl
This file contains the function called goldstein that implements the goldstein line search. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- ∇f (function) the gradient of the objective function.
- d (vector) vector containing a descent direction from the current estimation.

Moreover, you can specify the initial steplenght the constant used on the evaluation of Goldstein conditions and maximum number of iterations allowed.

## Example: 

```julia
  include("testfunctions.jl")   # see testfunctions.jl code for more details
  include("goldstein.jl")       # see goldsteinDY.jl code for more details 
  ndim = 100                    # dimension
  x = rand(ndim)
  d = -gradpen_I(x)
  t1, ierror1, fn1 = goldstein(pen_I, gradpen_I, x, d)
  t2, ierror2, fn2 = goldstein(pen_I, gradpen_I, x, d, α = 0.5, β = 1e-5, σ = 0.2, maxiter = 500)
```

This function returns the following informations:

- α (Float64) steplength which satisfies the goldstein conditions.
- ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
- fn (Int) number of function evaluations.

# golsdsteinDY.jl
This file contains the function called goldDY that implements the DY method (conjugate gradient) where the steplength holds the Goldstein conditions. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- gf (function) the gradient of the objective function.

Moreover, you can specify the tolerance the constant used on the evaluation process that will decide the direction taken and the maximum number of iterations allowed.

## Example: 

```julia
  include("testfunctions.jl")   # see testfunctions.jl code for more details
  include("goldsteinDY.jl")      # see cautiousDY.jl code for more details 
  ndim = 100                    # dimension
  x = rand(ndim)
  x1, fx1, normx1, iter1, ierror1, counter1, fn1, X1, Y1, Z1 = goldDY(x, pen_I, gradpen_I)
  x2, fx2, normx2, iter2, ierror2, counter2, fn2, X2, Y2, Z2 = goldDY(x, pen_I, gradpen_I, ϵ = 1.e-6, ϵ1 = 0.10, maxiter=500)
```

This function returns the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- fx (Float64) contain the value of objective function avaliated in the current estimation.
- normx (Float64) contain the value of the euclidean norm of the current estimation.
- iter (Int) number of iterations.
- ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
- counter (Int) number of the times where the gradient method was chosen.
- fn (Int) number of function evaluations.
- X (vector) contains the first coordinates in each estimation.
- Y (vector) contains the second coordinates in each estimation.
- Z (vector) contains the value of the function evaluation in each estimation.

## Remarks:

 - If error = 0, then a minimum was found.
 - If ndim = 2, you can use X, Y, Z to plot the sequences on level curves, for example.

# goldsteinPRP.jl
This file contains the function called goldsteinPRP that implements the PRP method (conjugate gradient) or the gradient method where the steplength holds the Goldstein conditions. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- ∇f (function) the gradient of the objective function.
- n (Int) dimension of the objective function.

Moreover, you can specify the maximum number of iterations allowed the tolerance and the method.

## Example: 

```julia
  include("testfunctions.jl")   # see testfunctions.jl code for more details
  include("goldsteinPRP.jl")       # see armijoPRP.jl code for more details
  ndim = 100                    # dimension
  x = rand(ndim)
  x1, fx1, normx1, iter1, ierror1, counter1, fn1, X1, Y1, Z1 = goldsteinPRP(x, pen_I, gradpen_I, ndim)
  x2, fx2, normx2, iter2, ierror2, counter2, fn2, X2, Y2, Z2 = goldsteinPRP(x, pen_I, gradpen_I, ndim, maxk = 400, ϵ = 1.e-6, method = 0)
```

This function returns the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- fx (Float64) contain the value of objective function avaliated in the current estimation.
- normx (Float64) contain the value of the euclidean norm of the current estimation.
- iter (Int) number of iterations.
- ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
- counter (Int) number of the times where the gradient method was chosen.
- fn (Int) number of function evaluations.
- X (vector) contains the first coordinates in each estimation.
- Y (vector) contains the second coordinates in each estimation.
- Z (vector) contains the value of the function evaluation in each estimation.

 ## Remarks:

 - If error = 0, then a minimum was found.
 - If ndim = 2, you can use X, Y, Z to plot the sequences on level curves, for example.
 - $maxk \in \mathbb{N}$.
 - $ϵ > 0$.
 - If method $\ne$ 1, then the gradient method will be choosed.

# performance_profile.jl
This file contains the function called performance_profile that plots the performance profile. To invoke this function you need the following informations:

- v (matrix) containing the the values of time (or any other metric) of a problem i and solver j.
- S (Vector{String}) contains the labels used.

Moreover you can specify the title.

## Example 1:

```julia
  include("performance_profile.jl")   # see performance_profile.jl code for more details
  x = rand(25, 3)                     # 25 problems and 3 solvers
  performance_profile(x, ["a"; "b"; "c"])
```
![alt text](https://github.com/wmd23/conjugategradient/blob/working/teste.png?raw=true)

## Example 1:

```julia
  include("performance_profile.jl")   # see performance_profile.jl code for more details
  x = rand(25, 3)                     # 25 problems and 3 solvers
  performance_profile(x, ["a"; "b"; "c"], t = "test_title")
```
![alt text](https://github.com/wmd23/conjugategradient/blob/working/teste1.png?raw=true)

## Remarks

- the number of solvers and the number of elements in the vector S must match (because each solver will have its label).

# PRPmethod.jl
This file contains the function called conjugadoPRP that implements the PRP method (conjugate gradient) where the steplength is computated by strong wolfe conditions. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- ∇f (function) the gradient of the objective function.
- n (Int) dimension of the objective function.

Moreover, you can specify the tolerance and the maximum number of iterations allowed.

## Example: 

```julia
  include("testfunctions.jl")   # see testfunctions.jl code for more details
  include("PRPmethod,jl")       # see PRPmethod.jl code for more details
  ndim = 100                    # dimension
  x = rand(ndim)
  x1, fx1, normx1, iter1, ierror1, counter1, fn1, X1, Y1, Z1 = conjugadoPRP(x, pen_I, gradpen_I, ndim)
  x2, fx2, normx2, iter2, ierror2, counter2, fn2, X2, Y2, Z2 = conjugadoPRP(x, pen_I, gradpen_I, ndim, maxk = 400, ϵ = 1.e-6, method = 0)
```

This function returns the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- fx (Float64) contain the value of objective function avaliated in the current estimation.
- normx (Float64) contain the value of the euclidean norm of the current estimation.
- iter (Int) number of iterations.
- ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
- counter (Int) number of the times where the gradient method was chosen.
- fn (Int) number of function evaluations.
- X (vector) contains the first coordinates in each estimation.
- Y (vector) contains the second coordinates in each estimation.
- Z (vector) contains the value of the function evaluation in each estimation.

 ## Remarks:

 - If error = 0, then a minimum was found.
 - If ndim = 2, you can use X, Y, Z to plot the sequence on level curves, for example.
