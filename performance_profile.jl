#
# Reference: E. Dolan and J. Moré, Benchmarking Optimization Software with Performance Profiles,
# mathematical Programming 91, pages 201-213, 2002. DOI 10.1007/s101070100263.
#

using Plots; gr()                   # GR backend

# function performance_ratio(V)
#-------------------------------
# function performance_ratio calculates the performance ratio

# Input Parameters
# =================
# V (matrix) containing the the values of time (or any other metric) of a problem i and solver j.

# Variables
# ==========
# x (Int) number of rows of the the matrix V.
# v (Int) number of columns of the the matrix V.
# min_column (Float64) lesser value in a specified column.

# Output
# ======
# A (matrix) contains the performance ratio of a problem i and solver j.

function performance_ratio(V)       # return the performance ratio matrix
    x = size(V)[1];                 # number of rows
    y = size(V)[2];                 # number of columns

    A = zeros(x,y)                  # matrix with desirable dimensions         

    for j in 1:y
        min_column = minimum(V[1:x, j])
        for i in 1:x
            A[i, j] = V[i, j]/min_column    # performance ratio formula
        end
    end

    return A

end

# function size_less(V, j, k)
#-----------------------------------------------
# function size_less calculates how many elements in a specific column is lesser than or equal to k

# Input Parameters
# =================
# V (matrix) containing the the values of performance ratio.
# j (Int) indicates the column.
# k (Float64) is used to make comparisons.

# Variables
# ==========
# counter (Int) number of elements in jth column that is lesser than or equal to k.

# Output
# ======
# counter (Int) number of elements in jth column that is lesser than or equal to k.

function size_less(V, j, k) 
    counter = 0

    for i in 1:size(V)[1]
        if V[i,j] ≤ k
            counter += 1
        end
    end

    return counter

end

function biggest(V)             # Calculate the biggest value in a matrix, except the value Inf
    big = 0
    for i in 1:size(V)[1]
        for j in 1:size(V)[2]
            if V[i,j] != Inf
                if V[i, j] > big
                    big = V[i, j]
                end
            end
        end
    end

    return big

end

function power_ten(k)       # Calculate the lesser power of ten which is greater than or equal to k
    t = 10
    while t < k
        t *= 10
    end
    
    return t

end

# function probability(V, s, t)
#-----------------------------------------------
# function probability implements the function which is plotted in a Performance Profile

# Input Parameters
# =================
# V (matrix) containing the the values of time (or any other metric) of a problem i and solver j.
# s (Int) indicates the column.
# t (Float64) is used to make comparisons.

# Variables
# ==========
# X (matrix) contains the performance ratio of a problem i and solver j.
# np (Int) np means number of problems.
# n (Inter) number of elements which is lesser than or equal to t.

# Output
# ======
# n/np (Float64) is a value between 0 and 1, see the reference for more details.

function probability(V, s, t)
    X = performance_ratio(V)
    np = size(X)[1]
    n = size_less(X, s, t)
    return n/np
end

# function performance_profile(v, S; t="")
#-----------------------------------------------
# function performance_profile plots the Performance Profile

# Input Parameters
# =================
# v (matrix) containing the the values of time (or any other metric) of a problem i and solver j.
# S (Vector{String}) contains the labels used.

# Optional Input Parameters
# =========================
# t (String) title of Performance Profile. (default = "")

# Variables
# ==========
# a (matrix{Float64}) contains the performance ratio of a problem i and solver j.
# b (Float64) biggest value of a, desconsidering Inf.
# n (Inter) number of elements which is lesser than or equal to t.
# x (UnitRange{Int64}) axis x range.
# p (Plots.Plot{Plots.GRBackend}) Performance Profile plot.
# f_array (Array{function}) contains the probability function for each column.

# Output
# ======
# p (Plots.Plot{Plots.GRBackend}) Performance Profile plot.

function performance_profile(v, S; t="")
    a = performance_ratio(v)
    b = biggest(a)
    x = range(1, power_ten(b), length = 100)
    p = 0;

    f_array = Array{Function}(undef, size(v)[2])

    # adding up functions to the f_array
    for j in 1:size(v)[2]
        f(x) = probability(v, j, x)
        f_array[j] = f
    end

    for j in 1:size(v)[2]
        if j == 1
            plot(x, f_array[j], title = t, leg = :outerright, dpi=200, ylabel = "Percentual de Problemas Resolvidos", xlabel ="log(performance ratio)" ,xlims = (1, Inf), label = S[j], xscale = :log10, linewidth = 2, gridstyle = :dash, xaxis =(formatter=x->string(round(log10(x); digits = 2))), ytick = ([0.2, 0.4, 0.6, 0.8, 1.0], ["20%", "40%", "60%", "80%", "100%"]))
        else
            p = plot!(x, f_array[j], label = S[j], linewidth = 2)
        end
    end

    return p

end