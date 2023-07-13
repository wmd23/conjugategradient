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
  conjugadoPRP(x, pen_I, gradpen_I, ndim)
```

This function returns the following informations:

- x (vector) vector containing the current estimation to be minimizer.
- fx (Float64) contain the value of objective function avaliated in the current estimation.
- normx (Float64) contain the value of the euclidean norm of the current estimation.
- k (Int) number of iterations.
- ierror (Int) the value stored in this variable tells the following messages: 0 - OK!, 1 - the maximum number of iterations has been exceeded.
- counter (Int) number of the times where the gradient method was chosen.
- fn (Int) number of function evaluations.
- X (vector) contains the first coordinates in each estimation.
- Y (vector) contains the second coordinates in each estimation.
- Z (vector) contains the value of the function evaluation in each estimation.
