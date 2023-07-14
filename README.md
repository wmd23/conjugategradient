# Armijo.jl
This file contain the function called armijo that calcutes a steplenght that Armijo's rule holds. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- gradfx (Float64) the value of the gradient avaliated in the current estimation.
- d (vector) vector containing a descent direction from the current estimation.

Moreover, you can specify the inicial steplenght $t$, the contant used on the backtracking process and the constant $c_1$ used on the evaluation of Armijo's rule $f(x+t\cdot d)\le f(x) + t\cdot c1\cdot \nabla f(x)^T\cdot d$.

## Example:

```julia
  include("testfunctions.jl")   # see testfunctions.jl for more details
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

- As you might know, you can use as many optional parameters to invoke you function as you want.
- $t > 0$.
- $y, c1 \in (0, 1)$.

# PRPmethod.jl
This file contain the function called conjugadoPRP that implements the PRP method (conjugate gradient) where the steplength is computated by strong wolfe conditions. To invoke this function you need the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- f (function) objective function.
- ∇f (function) the gradient of the objective function.
- n (Int) dimension of the objective function.

Moreover, you can specify the tolerance and the maximum number of iterations allowed.

## Example: 

```julia
  include("testfunctions.jl")   # see testfunctions.jl for more details
  ndim = 100                    # dimension
  x = rand(ndim)
  x, fx, normx, iter, ierror, counter, fn, X, Y, Z = conjugadoPRP(x, pen_I, gradpen_I, ndim)
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
