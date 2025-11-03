"""
   f, Df, Ψ, ijkl = makeODEmodel(α, γ, J,K,L, g::Vector{Symmetry})

make an ODE model dx/dt = f(x,R) of plane Couette flow for symmetry subgroup
<g[1], g[2], ... > discretization limits -J ≤ |j| ≤ J, -K ≤ k ≤ K, and 0 ≤ l ≤ L,
where j,k are indices of Fourier modes (e.g. cos αjx sin γkz) and l is the index
for the wall-normal polynomials (l=0 is lowest order, l=L highest order).

return functions f(x) and Df(x) (matrix of partial derivatives), basis set Ψ,
and ijkl indices for the basis set.
"""
function makeODEModel(α, γ, J,K,L, g::Vector{Symmetry})
    
    # make a complete basis up to J,K,L discretization level
    ijklfull = basisIndexMap(J,K,L)
    Ψfull = makeBasisSet(α, γ, ijklfull, normalize=false);
    
    # select the subset of those that are symmetric with all generators g[n]
    nsymm = findall([symmetric(ijklfull[n,:], g) for n in 1:size(ijklfull,1)])
    Ψ = Ψfull[nsymm]
    ijkl = ijklfull[nsymm,:]

    Nmodes = length(Ψ)
    y = Polynomial([0.0;1.0], :y)

    # these matrix calculations could be done with fewer conversions
    # should make elements with Rational{T} or Float type as desired
    # and matrices should be sprase from the start
    println("Making matrices B,A1,A2,S3...")
    B = [innerproduct(Ψ[i], Ψ[j]) for i in 1:Nmodes, j in 1:Nmodes]
    A1 = [-innerproduct(Ψ[i], y*xderivative(Ψ[j])) for i in 1:Nmodes, j in 1:Nmodes]
    A2 = [-innerproduct(Ψ[i], vex(Ψ[j])) for i in 1:Nmodes, j in 1:Nmodes]
    A3 = [innerproduct(Ψ[i], laplacian(Ψ[j])) for i in 1:Nmodes, j in 1:Nmodes]

    println("Making quadratic operator N...")
    Ndense = fill(0.0, Nmodes, Nmodes, Nmodes)
    for j in 1:Nmodes
        print("$j ")
        for k in 1:Nmodes    
            Ψj_dotgrad_Ψk = dotgrad(Ψ[j], Ψ[k])
            for i in 1:Nmodes
                val = -innerproduct(Ψ[i], Ψj_dotgrad_Ψk)
                Ndense[i,j,k] = abs(val) > 1e-15 ? val : 0
            end
        end
    end
    println("")
    flush(stdout)

    # convert rational matrices and dense N to floats and sparse float N
    Bf   = Float64.(B)
    A12f = Float64.(A1+A2)
    A3f  = Float64.(A3)
    N  = SparseBilinear(Ndense)
    Nf = SparseBilinear(N.ijk, Float64.(N.val), Nmodes)

    # precompute LU decomp of Bf to speed repeated calls to Bf x = b solves
    BfLU = lu(Bf)
    
    # construct Re-parameterized f and Df functions 
    f(x,R)  = BfLU\(A12f*x + (1/R)*(A3f*x) + Nf(x))
    Df(x,R) = BfLU\(A12f + (1/R)*(A3f) + derivative(Nf,x))
    
    f, Df, Ψ, ijkl # return functions f(x,R) and Df(x,R)
end

makeODEModel(α, γ, J,K,L) = makeODEmodel(α,γ,J,K,L, fill(Symmetry(), 0))
