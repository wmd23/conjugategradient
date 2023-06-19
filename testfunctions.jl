#
# Functions to be minimized and its gradients
#

using LinearAlgebra

# ackley function
function ackley(x, a=20, b=0.2, c=2π)
    d = length(x)
    return -a*exp(-b*sqrt(sum(x.^2)/d))-exp(sum(cos.(c*xi) for xi in x)/d) + a + exp(1)
end

# booth function
booth(x) = (x[1]+2x[2]-7)^2 + (2x[1]+x[2]-5)^2
gradbooth(x) = [2(x[1]+2x[2]-7)+4(2x[1]+x[2]-5); 4(x[1]+2x[2]-7) + 2(2x[1]+x[2]-5)]

#branin function
function branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π))
    return a*(x[2] - b*x[1]^2 + c*x[1] - r)^2 + s*(1 - t)*cos(x[1]) +s
end


#circle function
function circle(x)
    θ = x[1]
    r = 0.5 + 0.5*(2x[2]/(1+x[2]^2))
    y1 = 1 - r*cos(θ)
    y2 = 1 - r*sin(θ)
end

# flower function
function flower(x; a=1, b=1, c=4)
    return a*norm(x) + b*sin(c*atan(x[2], x[1]))
end

# michalewics function
function michalewics(x; m =10)
    return -sum(sin(v)*sin(i*v^2/π)^(2m) for (i, v) in enumerate(x))
end

# rosembrock function
rosembrock(x; a=1, b=2) = (a-x[1])^2 + b*(x[2] - x[1]^2)^2
gradrosembrock(x; a=1, b=2)=[-2*(a-x[1])-4*x[1]*b*(x[2] - x[1]^2); 2*b*(x[2] - x[1]^2)]

# wheeler's ridge function
wheeler(x, a=1.5) = -exp(-(x[1]*x[2] - a)^2 - (x[2] - a)^2)