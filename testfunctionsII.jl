using LinearAlgebra

# Rosenbrock function OK!
function rosenbrock(x)
    return 100*(x[2]-x[1]^2)^2 + (1-x[1])^2
end

function gradrosenbrock(x)
    g1 = -400*x[1]*(x[2]-x[1]^2) - 2*(1-x[1])
    g2 = 200*(x[2]-x[1]^2)
    return [g1; g2]
end

# Freudenstein and Roth function OK!
function freudenstein(x)
    return (-13+x[1]+((5-x[2])*x[2]-2)*x[2])^2 + (-29+x[1]+((x[2]+1)*x[2]-14)*x[2])^2
end

function gradfreudenstein(x)
    g1 = 2*(-13+x[1]+((5-x[2])*x[2]-2)*x[2]) + 2*(-29+x[1]+((x[2]+1)*x[2]-14)*x[2])
    g2 = 2*(-13+x[1]+((5-x[2])*x[2]-2)*x[2])*(-3*x[2]^2+10*x[2]-2) + 2*(-29+x[1]+((x[2]+1)*x[2]-14)*x[2])*(3*x[2]^2+2*x[2]-14)
    return [g1; g2]
end

# Powell badly scaled function (not OK!)

function powell(x)
    return (1e4*x[1]*x[2]-1)^2 + (exp(-x[1])+exp(-x[2])-1.0001)^2
end

function gradpowell(x)
    g1 = 2e4*x[2]*(1e4*x[1]*x[2]-1) - 2*exp(-x[1])*(exp(-x[1])+exp(-x[2])-1.0001)
    g2 = 2e4*x[1]*(1e4*x[1]*x[2]-1) - 2*exp(-x[2])*(exp(-x[1])+exp(-x[2])-1.0001)
    return [g1; g2]
end

# Brown badly scaled function (OK!)
function brown(x)
    return (x[1]-10^6)^2 + (x[2]-2*10^-6)^2+(x[1]*x[2]-2)^2
end
function gradbrown(x)
    g1= 2*(x[1]-10^6)+2*x[2]*(x[1]*x[2]-2)
    g2= 2*(x[2]-2*10^-6)+2*x[1]*(x[1]*x[2]-2)
    return  [g1;g2]
end

# Beale function (OK!)
function beale(x)
    return (1.5 - x[1]*(1-x[2]))^2 + (2.25 - x[1]*(1-x[2]^2))^2 + (2.625 - x[1]*(1-x[2]^3))^2
end

function gradbeale(x)
    g1 = -2*(1-x[2])*(1.5 - x[1]*(1-x[2])) - 2*(2.25 - x[1]*(1-x[2]^2))*((1-x[2]^2)) - 2*(2.625 - x[1]*(1-x[2]^3))*(1-x[2]^3)
    g2 = 2*x[1]*(1-x[2])*(1.5 - x[1]*(1-x[2])) + 4*x[1]*x[2]*(2.25 - x[1]*(1-x[2]^2)) + 6*x[1]*x[2]^2*(2.625 - x[1]*(1-x[2]^3))
    return [g1; g2] 
end

# Jennrich and Sampson function (not so sure!)
function jerrinch(x)
    return (2+2*im-exp(im*x[1])-exp(im*x[2]))^2
end

function gradjerrinch(x)
    g1 = 2*(2+2*im-exp(im*x[1])-exp(im*x[2]))*(-im*exp(im*x[1]))
    g2 = 2*(2+2*im-exp(im*x[1])-exp(im*x[2]))*(-im*exp(im*x[2]))
    return [g1;g2]
end

# Helical valley function (OK!)
function θ(x1, x2)
    if x1 > 0
        return 1/(2*pi)*atan(x2/x1)
    elseif x1 < 0
        return 1/(2*pi)*atan(x2/x1) + 0.5
    end
end

function helical(x)
    return 100*(x[3]-10*θ(x[1], x[2]))^2 + 100(sqrt(x[1]^2+x[2]^2)-1)^2 + x[3]^2
end

function gradhelical(x)
    g1 = 200*(sqrt(x[1]^2+x[2]^2)-1)*(x[1]/sqrt(x[1]^2+x[2]^2))
    g2 = 200*(sqrt(x[1]^2+x[2]^2)-1)*(x[2]/sqrt(x[1]^2+x[2]^2))
    g3 = 200*(x[3]-10*θ(x[1], x[2])) + 2*x[3]^2
    return [g1; g2; g3]
end

# Variably dimensioned function
function vardim(x)
    r = 0;
    s = 0;
    t = 0;
    for i in 1:length(x)
        r = (x[i] - 1)^2
        s = s + i*(x[i]-1)
        t = t + i*(x[i]-1)
    end

    s = s^2;
    t = t^4;

    return r+s+t
end

push!

function gradvardim(x)
    v = Float64[];
    s = 0;
    for i in 1:length(x)
        s = s + i*(x[i]-1)
    end

    for i in 1:length(x)
        push!(v, 2*(x[i]-1)+2*s*(x[i]-1+i)+4*s^3*(x[i]-1+i))
    end

    return v
end

# Extended Powell singular function

function singx(x)
    r = 0
    for i in 1:length(x)
        if 4%i == 1
            r = r + (x[i] + 10*x[i+1])^2
        elseif 4%i == 2
            r = r + 25*(x[i+1]-x[i+2])
        elseif 4%i == 3
            r = r + (x[i-1]-2*x[i])^4
        else
            r = r + 10*(x[i-3]-x[i])^4
        end
    end
end

function gradsingx(x)
    g = Float64[]
    for i in 1:lenght(x)
        if 4%i == 1
            push!(v, 2*(x[i]+10*x[i+1]) + 40*(x[i]-x[i+3])^3)
        elseif 4%i == 2
            push!(v, 20*(x[i-1]+10*x[i])+4*(x[i]-2*x[i+1])^3)
        elseif 4%i == 3
            push!(v, 10*(x[i]-x[i+1])-8*(x[i-1]-2*x[i]))
        else
            push!(v, -10*(x[i-1]-x[i])-40*(x[i-3]-x[i]))
        end
    end
end