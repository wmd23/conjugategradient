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
