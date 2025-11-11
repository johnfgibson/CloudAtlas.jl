function printflush(string)
    println(string)
    flush(stdout)
end

Base.@kwdef struct SearchParams{T<:Real}
    R::T=250.0
    ftol::T=1e-08 
    xtol::T=1e-08 
    δ::T=0.02
    Nnewton::Int=20 
    Nhook::Int=4
    Nmusearch::Int=6
    verbosity::Int=0
end

"""
hookstep(fx, Dfx, Δx, δ; Nmusearch=10, verbose=true, δtol = 0.01) :

    return hookstep Δx that minimizes 1/2 ||fx + Df Δx||^2 subject to ||Δx|| = δ
"""
function hookstep(fx, Dfx, δ, Δx_newt; δtol = 1e-04, Nmusearch=10, verbosity=0)
    
    norm_Δx = norm(Δx_newt)
    if (norm_Δx <= δ)
        verbosity > 0 && printflush("Newton step is within trust region, so returning hookstep Δx == Newton step")
        return Δx_newt
    end

    # start Newton iteration to solve norm(Δx(μ)) - δ == 0, with guess μ=0, Δx(μ) = Δx_newt
    μ = 0
    Δx = copy(Δx_newt)
    H = Dfx'*Dfx 
    
    verbosity > 0 && printflush("Starting search for hookstep Δx(μ) of radius norm(Δx(μ)) = δ = $δ")
    verbosity > 0 && printflush("μ = $(μ), norm(Δx(μ)) = $(norm_Δx)")
    
    # calculate hookstep Δx of radius δ using Newton iteration over μ on
    # equation norm(Δx(μ))- δ == 0. Read the algorithm notes above.
    
    for m=1:Nmusearch 
        verbosity > 1 && printflush("\nhookstep μ search $m")
        verbosity > 2 && printflush("Δx = $Δx")  
        ϕ = norm_Δx - δ           
        
        Hsolve = (H + μ*I)\Δx
        verbosity > 2 && printflush("H = $(H)")
        verbosity > 1 && printflush("μ = $μ")
        verbosity > 2 && printflush("H + μ I = $(H + μ*I)")
        verbosity > 2 && printflush("Hsolve = $Hsolve")
           
        ϕ = norm_Δx - δ           
        #ϕprime = -(1/norm_Δx) * Δx' * ((H + μ*I)\Δx)
        ϕprime = -(1/norm_Δx) * dot(Δx, (H + μ*I)\Δx)
        verbosity > 2 && printflush("ϕ  = $ϕ")
        verbosity > 2 && printflush("ϕ' = $ϕprime")               
        
        μ = μ - (norm_Δx/δ)*(ϕ/ϕprime)
        verbosity > 1 && printflush("new μ = $μ")
        
        Δx = -(H + μ*I)\(Dfx'*fx)
        
        verbosity > 2 && printflush("Dfx'*fx = $(Dfx'*fx)")
        verbosity > 2 && printflush("Δx = $Δx")
        norm_Δx = norm(Δx)
        verbosity > 1 && printflush("μ = $(μ), norm(Δx(μ)) = $(norm_Δx)")
                     
        if abs(norm_Δx - δ)/δ <= δtol
            verbosity > 0 && printflush("Found hookstep Δx with norm(Δx) = $(norm_Δx) ≈ $(δ) = δ.")
            return Δx
        end
    end
    verbosity > 0 && printflush("Stopping with hookstep Δx with norm(Δx) = $(norm_Δx) for $(δ) = δ.")
    return Δx
end


function Df_finitediff(f,x; eps=1e-06)
    fx = f(x)
    M = length(fx)
    N = length(x)
    
    Df = zeros(M,N)
    for j=1:N
        dxj = zeros(N)
        dxj[j] = eps
        fx_dxj = f(x + dxj)
        dfdxj = (fx_dxj - fx)/eps
        for i=1:M
            Df[i,j] = dfdxj[i]
        end
    end
    Df
end

"""
    hookstepsolve(model::ODEModel, xguess; δ=0.1, Nnewton=20, Nhook=4, Nmusearch=6, verbosity=0) :

Solve f(x) = 0 for x using Newton-hookstep algorithm using finite-difference estimate of Df
f is a function f(x)
xguess is the initial guess for the solution

optional named arguments:
  δ is the initial trust-region radius
  Nnewton is the maximum number of Newton steps
  Nhook is the maximum number of hooksteps per Newton step
  Nmusearch is the maximum number of iterations to find a hookstep with length equal to δ
  verbosity == 0,1,2 gives none, terse, and verbose diagnostic printouts

usage: 
  x, success = hookstep(f,Df, δ=0.2, Nnewton=8)   # example with two named arguments

return value x is the solution, success is a boolean flag indicating success (convergence) or failure (local minimum)
"""
function hookstepsolve(
    model::ODEModel,
    xguess::AbstractVector{T},
    params = SearchParams()
) where {T<:Real}
    return hookstepsolve(model.f, model.Df, xguess, params)
end


hookstepsolve(f, xguess, params = SearchParams()) = hookstepsolve(f, x -> Df_finitediff(f,x), xguess, params)

"""
hookstepsolve(f, Df, xguess; δ=0.1, Nnewton=20, Nhook=4, Nmusearch=6, verbosity=0) :

    Solve f(x) = 0 for x using Newton-hookstep algorithm with user-suppplied derivative Df
    f and Df are functions f(x), Df(x). Df supplies the derivative of f evaluated at x (Df = [df_i/dx_j]). 
    xguess is the initial guess for the solution
    δ is the initial trust-region radius
    Nnewton is the maximum number of Newton steps
    Nhook is the maximum number of hooksteps per Newton step
    Nmusearch is the maximum number of iterations to find a hookstep with length equal to δ
    verbosity == 0,1,2 gives none, terse, and verbose diagnostic printouts
"""
function hookstepsolve(
    f, 
    Df, 
    xguess::AbstractVector{T},
    params = SearchParams()
) where {T<:Real}

    δmax = 1
    δmin = 0.001
    
    #norm2(x) = x'*x
    norm2(x) = dot(x,x)
    
    # x, rx change once per newton step, are constant throughout hookstep calculations
    x = xguess           
    Xiterates = zeros(params.Nnewton+1, length(x))
    Xiterates[1,:] = x
    
    verbosity = params.verbosity
    R = params.R
    ftol = params.ftol
    xtol = params.xtol
    Nhook = params.Nhook
    Nmusearch = params.Nmusearch
    δ = params.δ

    #debugging = true

    for n = 1:params.Nnewton
        verbosity > 0 && printflush("\n\n\n===================================")
        verbosity > 0 && printflush("\nNewton step $n :")
        verbosity > 2 && printflush("x = $x")
        
        # start Newton step from best computed values over all previous computations
        verbosity > 0 && println("evaluating f(x)...")
        fx = f(x)
        rx = 1/2*norm2(fx)

        if sqrt(2rx) <= ftol
            verbosity > 0 && printflush("stopping and exiting search since norm(f(x)) = $(norm(fx)) ≤ $ftol = ftol.")
            #return x, Xiterates[1:n,:]
            return x, true
        else
            verbosity > 0 && printflush("continuing search since norm(f(x)) = $(sqrt(2rx)) > $ftol = ftol.")        
        end
        
        # compute Newton step Δx
        verbosity > 0 && printflush("computing Df(x)...")
        Dfx = Df(x)
        verbosity > 0 && printflush("solving Df(x) Δx = -f(x)...")
        Δx = -Dfx\fx
        #verbosity > 0 && printflush("mopping up Δx solve...")
        norm_Δx = norm(Δx)
        DfΔx_newt= Dfx*Δx

        verbosity > 2 && printflush("Δx newton = $Δx") ; flush(stdout) 
        #verbosity > 0 && printflush("norm(Δx newton) = $(norm(Δx))"); flush(stdout) 
        #verify that residual is decreasing in direction of Newton step (it oughta be!)
        #@show  fx
        #@show typeof(fx)
        #@show  fx'
        #@show typeof(fx')
        #@show (Dfx*Δx)
        #@show typeof(Dfx*Δx)
        #@show  fx'*(Dfx*Δx)
        #@show  typeof(fx'*(Dfx*Δx))
        
        #@show  dot(fx, Dfx*Δx)
        if  dot(fx, Dfx*Δx) >= 0 
            verbosity > 1 && printflush("Residual is increasing in direction of Newton step, indicating that the")
            verbosity > 1 && printflush("solution of the Newton-step equation is inaccurate. Exiting search and")
            verbosity > 1 && printflush("returning current value of x")
            Xiterates[n+1,:] = x + Δx
            #return x, Xiterates[1:n+1,:]
            return x, false
        end
             
        if norm(Δx) < xtol
            verbosity > 0 && printflush("Stopping because norm(Δx) = $(norm(Δx)) < $(xtol) = xtol")
            Xiterates[n+1, :] = x+Δx
            #return x+Δx, Xiterates[1:n+1,:]
            return x+Δx, false
        end
                
        Δx_newt = Δx           # Store Newton step and its norm for use in hookstep calculations       
        norm_Δx_newt = norm_Δx 
        x_hook = x + Δx        # Declare x_hook, modify it iteratively in hookstep loop.

        if verbosity > 0 
            relation = norm_Δx_newt ≤ δ ? "inside" : "outside"
            println("Starting trust region evaluation");
            println("Newton step radius  == $norm_Δx_newt")
            println("Trust region radius == $δ")
            println("Newton step is $relation the trust region");
        end

        for h in 1:Nhook

            if  norm_Δx_newt > δ
                verbosity > 0 && printflush("Hookstep $(h): finding hookstep Δx s.t. |Δx| ≤ δ = $δ")
                Δx = hookstep(fx, Dfx, δ, Δx_newt, Nmusearch=Nmusearch, verbosity=verbosity-1)
                verbosity > 1 && printflush("Found hookstep Δx s.t. |Δx| = δ = $δ")
            end

            verbosity > 1 && printflush("Computing Df Δx...")
            DfΔx = Dfx*Δx
            #δ = norm(Δx)  # hookstep function returns Δx with norm(Δx) ≈ δ, revise δ to make this exact
            
            if norm(Δx) < xtol
                verbosity > 0 && printflush("Stopping search because norm(Δx) = $(norm(Δx)) < $(xtol) = xtol")
                Xiterates[n+1, :] = x_hook
                #return x + Δx, Xiterates[1:n+1,:]
                return x + Δx, false
            end
                        
            newton_step_within_delta = norm_Δx_newt <= δ ? true : false

            verbosity > 0 && println("Assessing Δx...")
            # Compute actual (squared) residual of hookstep and linear & quadratic estimates based purely
            # on Δx and evaluations of f(x) and Df(x) at current Newton step. These derive from 
            # r(x + Δx) = 1/2 ||f(x+Δx)||^2
            #           ≈ 1/2 ||f(x) + Df(x) Δx||^2                  (this estimate is quadratic in Δx)
            #           ≈ 1/2 (f(x) + Df(x) Δx)ᵀ (f(x) + Df(x) Δx)
            #           ≈ 1/2 (f(x)ᵀ f(x) + 2 fᵀ Df(x) Δx + (Df(x) Δx)ᵀ (Df(x) Δx))
            #           ≈ r(x) + fᵀ Df(x) Δx + 1/2 (Df(x) Δx)ᵀ (Df(x) Δx) 
            #  
            # r(x + Δx) ≈ r(x) + fᵀ Df(x) Δx  (dropping quadratic terms give estimate linear in Δx)
            x_hook = x + Δx
            r_hook = 1/2*norm2(f(x+Δx))             # actual residual of hookstep, x + Δx
            #r_linear = 1/2*(norm2(fx) + 2*fx'*DfΔx)   # estimate of residual that is linear in Δx
            r_linear = rx + dot(fx,DfΔx)            # estimate of residual that is linear in Δx
            r_quadratic = 1/2*norm2(fx + DfΔx)      # estimate of residual that is quadratic in Δx

            #println("r(x)     = $rx")
            #println("r(x+Δxₕ) = $r_hook")
            #println("rlinear  = $r_linear")
            #println("rquadrat = $r_quadratic")

            # Differences in actual, linear, and quadratic estimates with sign set so that 
            # positive Δr == good, and bigger Δr == better. 
            Δr_hook = rx- r_hook
            #Δr_linear = (fx'*DfΔx)              # == -(r_linear - rx) without doing the subtraction
            Δr_linear = rx - r_linear            # dot(fx,DfΔx) is same without subtraction
            Δr_quadratic = rx - r_quadratic
   
            verbosity > 1 && printflush("\nTrust region adjustment for $(h): |Δx| = $δ, r(xₙ+Δx) = $(r_hook), compared to r(xₙ) = $(rx)")
            # revise trust region radius and do or don't recompute hookstep based on 
            # comparisons between actual change in residual and linear & quadratic models
            # note that diff between r_quadratic and r_linear is positive definite, so 0 < Δr_quadratic < Δr_linear
     
            #=====================================       
            if debugging
                hookdata = [rx; δ; r_linear; r_quadratic; r_hook]
                @show rx
                @show δ
                @show r_linear
                @show r_quadratic
                @show r_hook
                @show Δr_linear
                @show Δr_quadratic
                @show Δr_hook       
                m = -Δr_linear/δ
                c = (r_hook - rx - m*δ)/δ^2 
                δnew = -m/(2c)
                @show m
                @show c
                @show δnew
                save(hookdata, "hookdata")
                println("press enter to continue, x to stop debugging..")
                keystroke = readline()
                if keystroke == "x"
                    debugging = false
                end
            end
            ==============================#

            if Δr_hook > Δr_linear                # actual is better than linear estimate (quadratic helps!)
                verbosity > 0 && print("negative curvature, ")
                if newton_step_within_delta
                    verbosity > 0 && printflush("but newton step is within trust region")
                    verbosity > 0 && printflush("so don't increase δ, and go to next Newton step")
                    break
                else
                    verbosity > 0 && printflush("and newton step is outside trust region")
                    verbosity > 0 && printflush("so increase δ = $δ - > 3δ/2 = $(3δ/2) and recompute hookstep")
                    δ = 3δ/2           
                    continue
                end
            elseif Δr_hook < 0.01 * Δr_quadratic
                verbosity > 0 && printflush("poor improvement, decreasing δ  = $δ -> δ/2 = $(δ/2) and recomputing hookstep") 
                δ = δ/2                            # actual is vastly worse than quadratic estimate
                continue                           # reduce trust region  and recompute hookstep
            elseif Δr_hook < 0.10 * Δr_quadratic
                verbosity > 0 && printflush("marginal improvement, decreasing δ  = $δ -> δ/2 = $(δ/2) and continuing to next Newtown step")
                δ = δ/2                            # actual is worse than quadratic estimate
                break                              # reduce trust region and go to next Newton step
            elseif  0.8*Δr_quadratic < Δr_hook < 1.2*Δr_quadratic   
                # Actual decrease at current hookstep is close to quadratic estimate from Df,
                # and not as good as suggested by slope at origin (we have Δr_hook > Δr_linear, so r_hook >= r_linear).
                # Attempt to adjust trust-region radius δ based on a quadratic model of residual, r(δ)
                # using form r(δ) = rx + m δ + c δ^2 and the measured value of r(δ) = r_hook (where δ is the current radius)
                # Get m from conditions at δ=0, namely, m = fx'*DfΔx/||Δx|| = fx'*DfΔx/δ
                # Fit c to data point r(δ) = r_hook. That gives c = (r_hook - rx - m δ)/δ^2
                # We know c>0, since r_hook > r_linear = rx + + m δ. So the model has a minimum between
                # 0 and the current δ. The minimum occurs at δ = -m/(2c)
                m = -Δr_linear/δ
                c = (r_hook - rx - m*δ)/δ^2 
                δnew = -m/(2c)
                #@show m
                #@show c
                #@show δnew
                verbosity > 0 && printflush("accurate improvement, quadratic model suggests δ  = $δ -> δnew = $(δnew)")
                if newton_step_within_delta
                    verbosity > 0 && printflush("but Newton step is within trust region, so don't change δ, and go to next Newton step")
                elseif δnew < 2δ/3
                    verbosity > 0 && printflush("too much decrease, changing δ = $δ -> 2δ/3 = $(2δ/3) and continuing to next Newtown step")
                    δ = 2δ/3
                elseif δnew > 4δ/3    
                    verbosity > 0 && printflush("too much increase, changing δ = $δ -> 4δ/3 = $(4δ/3) and continuing to next Newtown step")
                    δ = 4δ/3
                else
                    verbosity > 0 && printflush("not too much change, changing δ = $δ -> δnew = $δnew and continuing to next Newtown step")
                    δ = δnew
                end

                break
            else 
                verbosity > 0 && printflush("good improvement, keeping δ = $(δ) and continuing to next Newtown step")
                break                              # hookstep is decent enough, don't adjust, go to next Newton               
            end
        end
        
        verbosity > 0 && printflush("Finished trust region evaluation, trust region radius δ = $δ...")
        verbosity > 0 && printflush("norm(Δx hook) = $(norm(x - x_hook))")
        verbosity > 2 && printflush("prior x == $x")
        x = x_hook
        verbosity > 2 && printflush("  new x == $x == x_hook")
        Xiterates[n+1,:] = x
    end

    
    verbosity > 0 && printflush("Stopping because we reached maximum # of Newton iterations")
    #return x, Xiterates
    return x, false
end
