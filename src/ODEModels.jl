"""
   f, Df = ODEmodel(Ψ)

Make an ODE model dx/dt = f(x,R) of plane Couette flow by Galerkin projection
onto basis set Ψ. 

Return functions f(x,R) and Df(x,R) (the matrix of partial derivatives of f). 
"""
function ODEModel(Ψ::AbstractVector{BasisFunction{T}}) where {T<:Real}
    
    Nmodes = length(Ψ)
    y = Polynomial{T, :y}([zero(T), one(T)])

    # these matrix calculations could be done with fewer conversions
    # should make elements with Rational{T} or Float type as desired
    # and matrices should be sprase from the start
    println("Making matrices B,A1,A2,S3...")
    B = [innerproduct(Ψ[i], Ψ[j]) for i in 1:Nmodes, j in 1:Nmodes]
    A1 = [-innerproduct(Ψ[i], y*xderivative(Ψ[j])) for i in 1:Nmodes, j in 1:Nmodes]
    A2 = [-innerproduct(Ψ[i], vex(Ψ[j])) for i in 1:Nmodes, j in 1:Nmodes]
    A3 = [innerproduct(Ψ[i], laplacian(Ψ[j])) for i in 1:Nmodes, j in 1:Nmodes]

    println("Making quadratic operator N...")
    Ndense = zeros(T, Nmodes, Nmodes, Nmodes)
    for j in 1:Nmodes
        print("$j ")
        for k in 1:Nmodes    
            Ψj_dotgrad_Ψk = dotgrad(Ψ[j], Ψ[k])
            for i in 1:Nmodes
                val = -innerproduct(Ψ[i], Ψj_dotgrad_Ψk)
                Ndense[i,j,k] = abs(val) > 1e-15 ? val : zero(T) 
            end
        end
    end

    A12 = A1+A2
    N  = SparseBilinear(Ndense)

    # precompute LU decomp of Bf to speed repeated calls to Bf x = b solves
    Bfact = lu(B)
    @show Bfact
    
    # construct Re-parameterized f and Df functions 
    f(x,R)  = Bfact\(A12*x + (1/R)*(A3*x) + N(x))
    Df(x,R) = Bfact\(A12 + (1/R)*(A3) + derivative(N,x))
    
    f, Df # return functions f(x,R) and Df(x,R)
end

