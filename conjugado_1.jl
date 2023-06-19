# Armijo's line search
# f - function to be optimized
# ∇f - gradient of function f
# x - guess point
# d - descent direction
# γ ∈ (0,1) e c1 ∈ (0,1/2)

using Printf

function armijo(f, ∇f, x, d; γ=0.9, c1=1e-4, maxiter=1000)
    t = 1;
    iter = 0;
    gradfx = ∇f(x); 
    k = gradfx'*d
    k = k[1]
    ierror = 0

    while (f(x+t*d) > f(x) + c1*t*k) && (iter < maxiter)
        t = γ*t
        iter = iter + 1
    end

    if iter < maxiter
        return t,ierror
    else
        return t,ierror
    end
end

#
# Goldstein's Linesearch (adapted from the function below)
#
function goldstein(f, ∇, x, d; α=1, β=1e-4, σ=0.1, maxiter = 1000)
    y0, g0, y_prev, α_prev = f(x), ∇(x)⋅d, NaN, 0
    αlo, αhi = NaN, NaN
    ierror, iter = 0, 0

    #bracket phase
    while true
        y = f(x + α*d)
        if y0 + (1-β)*α*g0 > y || (!isnan(y_prev) && y ≥ y_prev)
            αlo, αhi = α_prev, α
            break
        end
        g = ∇(x + α*d)⋅d
        if y ≤ y0 + β*α*g0
            return α,ierror
        elseif g ≥ 0
            αlo, αhi = α, α_prev
            break
        end
        y_prev, α_prev, α = y, α, 2α
    end

    # zoom phase
    ylo = f(x + αlo*d)
    while true
        α = (αlo + αhi)/2
        iter=iter+1
        if iter>maxiter
            ierror = 3
            return α, ierror
        end
        y = f(x + α*d)
        if y0 + (1-β)*α*g0 > y || y0 + (1-β)*α*g0 ≥ ylo
            αhi = α
        else
            g = ∇(x + α*d)⋅d
            if y ≤ y0 + β*α*g0
                return α, ierror

            elseif g*(αhi - αlo) ≥ 0
                αhi = αlo
            end
            αlo = α
        end
    end
end

#
# Strong Wolfe's Linesearch 
# code found on the book Algorithms for Optimization by Mykel J. Kochenderfer and Tim A. Wheeler
#

function wolfe(f, ∇, x, d; α=1, β=1e-4, σ=0.1, maxiter = 1000)
    y0, g0, y_prev, α_prev = f(x), ∇(x)⋅d, NaN, 0
    αlo, αhi = NaN, NaN
    ierror, iter = 0, 0

    #bracket phase
    while true
        y = f(x + α*d)
        if y > y0 + β*α*g0 || (!isnan(y_prev) && y ≥ y_prev)
            αlo, αhi = α_prev, α
            break
        end
        g = ∇(x + α*d)⋅d
        if abs(g) ≤ -σ*g0
            return α, ierror
        elseif g ≥ 0
            αlo, αhi = α, α_prev
            break
        end
        y_prev, α_prev, α = y, α, 2α
    end

    # zoom phase
    ylo = f(x + αlo*d)
    while true
        α = (αlo + αhi)/2
        iter = iter + 1
        if iter > maxiter
            ierror = 3
            return α, ierror 
        end
        y = f(x + α*d)
        if y > y0 + β*α*g0 || y ≥ ylo
            αhi = α
        else
            g = ∇(x + α*d)⋅d
            if abs(g) ≤ -σ*g0
                return α, ierror

            elseif g*(αhi - αlo) ≥ 0
                αhi = αlo
            end
            αlo = α
        end
    end
end



# n - dimension
# linesearch is used to select which linesearch to used
# method == 0 (gradient method), method != 0 (conjugated gradient method)
# if you do not want to print the table, please coment the lines: 

function conjugado(x, f, ∇f, n, linesearch, method; maxiter = 1000, ϵ = 1e-6)
    d = -∇f(x)
    gradfx = -d
    ierror = 0
    iter = 0
    t = NaN

    #println("   iter|     ||∇f(x)||      |     t⋅||dk||      ")
    #println("-----------------------------------------------")

    t0 = time()

    while norm(∇f(x)) > ϵ && iter < maxiter

        #newgradfx = ∇f(x)
        #normgradfx = norm(newgradfx,2)
        #normdk = norm(d,2)

        #@printf("%5d  |%18.10e  |%18.10e\n",iter,normgradfx,t * normdk)

        t,error = linesearch(f, ∇f, x, d)

        x = x + t*d
        newgradfx = ∇f(x)
        if method == 0
            βk = 0
        else
            if (iter+1)%n != 0 
                βk = (newgradfx'⋅newgradfx)/(gradfx'⋅gradfx) # Fletcher and Reeves
            else
                βk = 0
            end
        end
        d = -newgradfx + βk*d
        iter = iter + 1
        gradfx = newgradfx
    end

    #newgradfx = ∇f(x)
    #normgradfx = norm(newgradfx,2)
    #normdk = norm(d,2)

    #@printf("%5d  |%18.10e  |%18.10e\n",iter,normgradfx,t * normdk)

    if iter < maxiter
        et = time() - t0
        return x, iter, et, ierror
    else
        ierror = 1
        et = time() - t0
        return x, iter, et, ierror
    end
end