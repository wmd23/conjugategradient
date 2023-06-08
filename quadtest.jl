using LinearAlgebra
include("quadcgrad.jl")
n = 1
A = rand(n,n)
A = 0.5 * (A + A')
diag = ones(n)
Id = diagm(diag)
A = A + n * Id
b = rand(n,1)
c = rand(1,1)

function f(x)
    z =  0.5 * x'*A*x+b'*x+c
    return z[1]
end

gf(x) = A*x+b
x0 = rand(n,1)
maxiter = 1000
ϵ = 1.e-8

(sol,iter) = quadcgrad(f,gf,x0,ϵ,maxiter)
