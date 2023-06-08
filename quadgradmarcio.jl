function congrad(x,f,∇f,ϵ,maxiter,A) 
    
    ierror=0
    iter=0
    d0 = -∇f(x)
    k=0
    tk=0.0
    dk=0.0
    β = 0.0
    fx = f(x)
    gradfx = ∇f(x)

    α = NaN

    while true
        newgradfx = ∇f(x)
        normgradfx = norm(newgradfx,2)
        normdk = norm(dk,2)

        @printf("%5d  %18.10e  %18.10e\n",iter,normgradfx,α * normdk)

        if normgradfx < ϵ
            return(x)
        end

        iter += 1
        if iter > maxiter
            ierror = 1
            return(x)
        end

        if k == 0
            tk = -(newgradfx' * d0/ d0'*A*d0)
            x = x+tk*d0
            gradfx = ∇f(x)
            β = (d0' * A * gradfx/d0'* A * d0)
            dk = - gradfx + β*d0
            k += 1
        else
            tk = -(newgradfx' * dk/ dk'*A*dk)
            x = x+tk*dk
            gradfx = ∇f(x)
            β = (dk' * A * gradfx/dk'* A * dk)
            dk = - gradfx + β*dk
            k += 1
        end
    end
end