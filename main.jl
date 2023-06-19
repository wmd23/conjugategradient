using LinearAlgebra, DataFrames, Random, Printf, Plots, BenchmarkProfiles

include("testfunctions.jl")
include("conjugado_1.jl")

nguess=5 #number of guess(es)

ndim = [2] #vector which contains the dimension number (only use if the number of dimensions is variable)

T=Float64[] #create an empty Float64 vector with time
I = Float64[] #create an empty Float64 vector with iterates

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
        x = rand(rng,2,1)
        (x, iter, et, ierror) = conjugado(x,rosembrock, gradrosembrock, 2, linesearch, 1)

        if ierror > 0
            push!(T,Inf)
            push!(I, Inf)
        else
            push!(T, et)
            push!(I, iter)
        end  
        #println("$n $etime")
    end
end

h=nguess*size(ndim,1)
Z=[I[1:h] I[h+1:2h] I[2h+1:3h]] #Matriz com os tempos
performance_profile(PlotsBackend(), Z, ["MGCA", "MGCG","MGCW"], title="Deus é bom!",xlabel="CPU time ratio",ylabel="Percentual de problemas resolvidos", legend=:bottomright, lw=2)
savefig("performanceprofile3.png")
#Z=[I[1:h] I[h+1:2h] I[2h+1:3h]] #Matriz com as iteradas
#performance_profile(PlotsBackend(), Z, ["MGCA", "MGCG", "MGCW"],xlabel="CPU time ratio",ylabel="Percentual de problemas resolvidos", legend=:bottomright, lw=2)
#savefig("performanceprofile2.png")