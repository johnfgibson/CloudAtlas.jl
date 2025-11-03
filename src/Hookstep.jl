function printflush(string)
    println(string)
    flush(stdout)
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
hookstepsolve(f, Df, xguess; δ=0.1, Nnewton=20, Nhook=4, Nmusearch=6, verbosity=0) :

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

function hookstepsolve(f, xguess; ftol=1e-08, xtol=1e-08, δ=0.1, Nnewton=20, Nhook=4, Nmusearch=6, verbosity=0)
    hookstepsolve(f, x -> Df_finitediff(f,x), xguess, ftol=ftol, xtol=xtol, δ=δ, Nnewton=Nnewton,
                  Nhook=Nhook, Nmusearch=Nmusearch, verbosity=verbosity)
end

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
function hookstepsolve(f, Df, xguess; ftol=1e-08, xtol=1e-08, δ=0.1, Nnewton=20, Nhook=4, Nmusearch=6, verbosity=0)

    δmax = 1
    δmin = 0.001
    
    #norm2(x) = x'*x
    norm2(x) = dot(x,x)
    
    # x, rx change once per newton step, are constant throughout hookstep calculations
    x = xguess           
    Xiterates = zeros(Nnewton+1, length(x))
    Xiterates[1,:] = x
    
    for n = 1:Nnewton
        verbosity > 0 && printflush("\nNewton step $n :")
        verbosity > 2 && printflush("x = $x")
        
        # start Newton step from best computed values over all previous computations
        verbosity > 0 && print("evaluating f(x)...")
        fx = f(x)
        rx = 1/2*norm2(fx)
        verbosity > 0 && printflush("norm(fx) = $(sqrt(2rx))")

        if sqrt(2rx) <= ftol
            verbosity > 0 && printflush("stopping and exiting search since norm(f(x)) = $(norm(fx)) ≤ $ftol = ftol.")
            #return x, Xiterates[1:n,:]
            return x, true
        else
            verbosity > 0 && printflush("continuing search since norm(f(x)) = $(sqrt(2rx)) > $ftol = ftol.")        
        end
        
        # compute Newton step Δx
        verbosity > 1 && printflush("computing Df(x)...")
        Dfx = Df(x)
        verbosity > 1 && printflush("solving Df(x) Δx = -f(x)...")
        Δx = -Dfx\fx
        verbosity > 1 && printflush("mopping up Δx solve...")
        norm_Δx = norm(Δx)
        DfΔx_newt= Dfx*Δx

        verbosity > 2 && printflush("Δx newton = $Δx") ; flush(stdout) 
        verbosity > 0 && printflush("norm(Δx newton) = $(norm(Δx))"); flush(stdout) 
        # verify that residual is decreasing in direction of Newton step (it oughta be!)
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

        for h in 1:Nhook
            verbosity > 1 && printflush("\nHookstep $(h): finding hookstep Δx s.t. |Δx| = δ = $δ")
            Δx = hookstep(fx, Dfx, δ, Δx_newt, Nmusearch=Nmusearch, verbosity=verbosity-1)
            verbosity > 1 && printflush("Found hookstep Δx s.t. |Δx| = δ = $δ")

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

            verbosity > 1 && printflush("Assessing Δx...")
            # Compute actual (squared) residual of hookstep and linear & quadratic estimates based purely
            # on Δx and evaluations of f(x) and Df(x) at current Newton step. These derive from 
            # r(x + Δx) = 1/2 ||f(x+Δx)||^2
            #           ≈ 1/2 ||f(x) + Df(x) Δx||^2                  (this estimate in quadratic in Δx)
            #           ≈ 1/2 (f(x) + Df(x) Δx)ᵀ (f(x) + Df(x) Δx)
            #           ≈ 1/2 (f(x)ᵀ f(x) + 2 fᵀ Df(x) Δx + (Df(x) Δx)ᵀ (Df(x) Δx))
            #           ≈ r(x) + fᵀ Df(x) Δx + 1/2 (Df(x) Δx)ᵀ (Df(x) Δx) 
            #  
            # r(x + Δx) ≈ r(x) + fᵀ Df(x) Δx  (dropping quadratic terms give estimate linear in Δx)
            x_hook = x + Δx
            r_hook = 1/2*norm2(f(x+Δx))             # actual residual of hookstep, x + Δx
            #r_linear = 1/2*(norm2(fx) + fx'*DfΔx)   # estimate of residual that is linear in Δx
            r_linear = 1/2*(norm2(fx) + dot(fx,DfΔx))   # estimate of residual that is linear in Δx
            r_quadratic = 1/2*norm2(fx + DfΔx)      # estimate of residual that is quadratic in Δx
        
            # Differences in actual, linear, and quadratic estimates with sign set so that 
            # positive Δr == good, and bigger Δr == better. 
            Δr_hook = -(r_hook - rx)
            #Δr_linear = -1/2*(fx'*DfΔx)              # == -(r_linear - rx) without doing the subtraction
            Δr_linear = -1/2*dot(fx,DfΔx)             # == -(r_linear - rx) without doing the subtraction
            Δr_quadratic = -(r_quadratic - rx)
                        
            verbosity > 1 && printflush("\nTrust region adjustment for $(h): |Δx| = $δ, r(xₙ+Δx) = $(r_hook), compared to r(xₙ) = $(rx)")
            # revise trust region radius and do or don't recompute hookstep based on 
            # comparisons between actual change in residual and linear & quadratic models
            # note that diff between r_quadratic and r_linear is positive definite, so 0 < Δr_quadratic < Δr_linear
            
            if Δr_hook > Δr_linear                # actual is better than linear estimate (quadratic helps!)
                verbosity > 0 && print("negative curvature, ")
                if newton_step_within_delta
                    verbosity > 1 && printflush("but newton step is within trust region")
                    verbosity > 1 && printflush("so don't increase δ, and go to next Newton step")
                    break
                else
                    verbosity > 1 && printflush("and newton step is outside trust region")
                    verbosity > 1 && printflush("so increase δ = $δ - > 3δ/2 = $(3δ/2) and recompute hookstep")
                    δ = 3δ/2           
                    continue
                end
            elseif Δr_hook < 0.01 * Δr_quadratic
                verbosity > 1 && printflush("poor improvement, decreasing δ  = $δ -> δ/2 = $(δ/2) and recomputing hookstep") 
                δ = δ/2                            # actual is vastly worse than quadratic estimate
                continue                           # reduce trust region  and recompute hookstep
            elseif Δr_hook < 0.10 * Δr_quadratic
                verbosity > 1 && printflush("marginal improvement, decreasing δ  = $δ -> δ/2 = $(δ/2) and continuing to next Newtown step")
                δ = δ/2                            # actual is worse than quadratic estimate
                break                              # reduce trust region and go to next Newton step
            elseif 0.8*Δr_quadratic < Δr_hook < 2* Δr_quadratic   
                # actual is close to quadratic estimate                
                # revise trust region based on quadratic model of residual and recompute
                rprime = Δr_linear/δ
                rprime2 = 2*(r_hook - rx - rprime*δ)/δ^2 
                δnew = -rprime/rprime2
                verbosity > 1 && printflush("accurate improvement, considering δ  = $δ -> δnew = $(δnew) from quadratic model")
                                      
                if newton_step_within_delta
                    verbosity > 1 && printflush("but Newton step is within trust region, so don't change δ, and go to next Newton step")
                elseif δnew < 2δ/3
                    verbosity > 1 && printflush("too much decrease, changing δ = $δ -> 2δ/3 = $(2δ/3) and continuing to next Newtown step")
                    δ = 2δ/3
                elseif δnew > 4δ/3    
                    verbosity > 1 && printflush("too much increase, changing δ = $δ -> 4δ/3 = $(4δ/3) and continuing to next Newtown step")
                    δ = 4δ/3
                else
                    verbosity > 1 && printflush("not too much change, changing δ = $δ -> δnew = $δnew and continuing to next Newtown step")
                    δ = δnew
                end
                break
            else 
                verbosity > 1 && printflush("good improvement, keeping δ = $(δ) and continuing to next Newtown step")
                break                              # hookstep is decent enough, don't adjust, go to next Newton               
            end
                    
        end
        
        verbosity > 1 && printflush("Finished with hookstep computations. Resetting x to x_hook")
        verbosity > 2 && printflush("prior x == $x")
        x = x_hook
        verbosity > 2 && printflush("  new x == $x == x_hook")
        Xiterates[n+1,:] = x
    end

    
    verbosity > 0 && printflush("Stopping because we reached maximum # of Newton iterations")
    #return x, Xiterates
    return x, false
end
