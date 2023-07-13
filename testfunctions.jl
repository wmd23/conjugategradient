#
# Reference: Moré, J.J., Garbow, B.S. and Hillstrome, K.E., 1981, Testing unconstrained optimization software.
# ACM Transactions on Mathematical Software, 7, 17-41.
#

# Penalty function I
function pen_I(x; a=1.e-5)
    r, s = 0, 0
    for i in 1:length(x)
        r += (x[i]-1)^2
        s += x[i]^2
    end
    s = (s-0.25)^2
    r = r*a
    return r + s
end

function gradpen_I(x; a=1.e-5)
    v = Float64[]
    s = 0
    for j in 1:length(x)
        s += x[j]^2
    end
    for i in 1:length(x)
        push!(v, 2*a*(x[i]-1)+4*x[i]*(s-0.25))
    end
    return(v)
end

#
# Reference: Surjanovic, S. & Bingham, D. (2013). Virtual Library of Simulation Experiments: 
# Test Functions and Datasets. Retrieved June 29, 2023, from http://www.sfu.ca/~ssurjano.
#

# Trid Function (Bowl-Shaped)

function trid(x)
    r, s = 0, 0
    for j in 1:length(x)
        s += (x[j]-1)^2

        if j != 1
            r += x[j]*x[j-1]
        end
    end

    return s - r

end

function gradtrid(x)
    v = Float64[]
    l = length(x)
    for i in 1:l
        if i == 1
            push!(v, 2*(x[i]-1) - x[i+1])
        elseif i == l
            push!(v, 2*(x[i]-1) - x[i-1])
        else
            push!(v, 2*(x[i]-1) - x[i-1] - x[i+1])
        end
    end

    return v

end