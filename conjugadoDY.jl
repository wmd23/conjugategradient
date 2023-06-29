function armijo(f, ∇f, x, d; γ=0.5, c1=0.1, maxk=1000)
    t = 1;
    k = 0;
    gradfx = ∇f(x); 
    k = gradfx'*d
    k = k[1]
    ierror = 0

    while (f(x+t*d) > f(x) + c1*t*k) && (k < maxk)
        t = γ*t
        k = k + 1
    end

    if k < maxk
        return t, ierror
    else
        ierror = 2
        return t, ierror
    end
end

function wolfe(f, ∇, x, d; α=1, β=1e-4, σ=0.1, maxk = 1000)
    y0, g0, y_prev, α_prev = f(x), ∇(x)⋅d, NaN, 0
    αlo, αhi = NaN, NaN
    ierror, k = 0, 0

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
        k = k + 1
        if k > maxk
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

function conjugadoPRP(x, f, ∇f, n; maxk = 10000, ϵ = 1e-5)
    d = -∇f(x)
    gradfx = -d
    ierror = 0
    k = 0
    counter = 1;

    while norm(∇f(x)) > ϵ && k < maxk

        t,error = wolfe(f, ∇f, x, d)
        if ierror > 0
            println(error)
        end
        x = x + t*d
        newgradfx = ∇f(x)
        if (k+1)%n != 0 
            βk = (newgradfx'⋅(newgradfx - gradfx))/(gradfx'⋅gradfx)
            βk = max(0, βk)
            if βk == 0
                counter = counter + 1
            end
        else
            βk = 0
            counter = counter + 1
        end
        d = -newgradfx + βk*d
        k = k + 1
        gradfx = newgradfx
    end

    if k < maxk
        fx = f(x)
        return x, fx, k, counter
    else
        ierror = 1
        fx = f(x)
        return x, fx, k, counter
    end
end

function conjugadoDY(x, f, ∇f; maxk = 10000, ϵ = 1e-5, ϵ1=0.01)
    d = -∇f(x)
    gradfx = -d
    newgradfx = gradfx
    ierror = 0
    k = 0;
    counter = 1;

    while norm(∇f(x)) > ϵ && k < maxk

        if k == 0
            d = -∇f(x)
        else
            if d⋅(newgradfx - gradfx) < ϵ1*norm(d)*norm(gradfx)
                gradfx = newgradfx
                d = -newgradfx
                counter = counter + 1
            else
                βk = norm(newgradfx)^2/d⋅(newgradfx - gradfx)
                d = -newgradfx + βk*d
                gradfx = newgradfx
            end
        end
        t, error = armijo(f, ∇f, x, d)
        # if ierror > 0
        #     println(error)
        # end
        x = x + t*d
        newgradfx = ∇f(x)
        k = k + 1
    end

    if k < maxk
        fx = f(x)
        return x, fx, k, counter
    else
        ierror = 1
        fx = f(x)
        return x, fx, k, counter
    end
end