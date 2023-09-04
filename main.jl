using LinearAlgebra, DataFrames, Random, Printf, Plots, BenchmarkProfiles

include("testfunctions.jl")
include("conjugado.jl")

nguess=5 #number of guess(es)

ndim = [2] #vector which contains the dimension number (only use if the number of dimensions is variable)

T=Float64[] #create an empty Float64 vector

rng = MersenneTwister(12345)

for B in 1:3 
    if B==1 
        println("Armijo")
        linesearch = armijo
    elseif B==2 
        println("Goldstein")
        linesearch = goldstein
    else
        println("Wolfe")
        linesearch = wolfe
    end
    for m in 1:nguess
        x = rand(snd,n,1)
        (x,ierror,info,etime) = conjugado(x, booth, gradbooth, 2, linesearch, 1)

        if ierror > 0
            push!(T,Inf)
        else 
            push!(T, etime)
        end   

        println("$n $etime")
    end
end

h=nguess*size(ndim,1)
Z=[T[1:h] T[h+1:2h] T[2h+1:3h]] #Matriz com os tempos
performance_profile(PlotsBackend(), Z, ["MGCA", "MGCG", "MGCW"],xlabel="CPU time ratio",ylabel="Percentual de problemas resolvidos", legend=:bottomright, lw=2)
savefig("performanceprofile.png")