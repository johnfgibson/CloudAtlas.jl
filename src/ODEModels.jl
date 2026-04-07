import Base.length

"""
ODEModel
"""
struct ODEModel{T<:Real, TB, TF, TDF}
    α::T                         # streamwise wavenumber α = 2pi/Lx
    γ::T                         # spanwise wavenumber γ = 2pi/Lz
    H::Vector{Symmetry}          # symmetry subgroup
    ijkl::Matrix{Int}            # set of allowed i,j,k,l indices 
    Ψ::Vector{BasisFunction{T}}  # basis set built from allowed i,j,k,l
    Ψshear::Vector{T}            # value of dΨudy at walls, used to calculate shear rate
    B::AbstractMatrix{T}         # inner product matrix 
    A1::AbstractMatrix{T}        # linear term arising from base flow
    A2::AbstractMatrix{T}        # linear term arising from Laplacian
    N::SparseBilinear{T}         # nonlinear term from u dot grad u
    Bfact::TB                    # LU factorization of inner product matrix
    f::TF                        # RHS function for ODE dx/dt = f(x,R) 
    Df::TDF                      # Df(x,R), the derivative of f, [Df]_ij = df_i/dx_j
end

function ODEModel(α::T, γ::T, J::Int, K::Int, L::Int, H::Vector{Symmetry}; normalize=false) where T<:Real
    ijkl = basisIndices(J,K,L,H)    # Compute set of allowed ijkl for H-symmetric Ψijkl 
    m = size(ijkl,1)                # State space dimension
    println("J,K,L,m == $J,$K,$L,$m")
    println("(2J+1)(2K+1)(2L+1) + 1 == $((2J+1)*(2K+1)*(2L+1) + 1)")
        
    Ψ = basisSet(α, γ, ijkl, normalize=normalize)  # Construct basis set 
    Ψshear = zeros(T, m)
    for i=1:m
        if Ψ[i].u[1].ejx.waveindex == 0 && Ψ[i].u[1].ekz.waveindex == 0
            dΨudy = yderivative(Ψ[i].u[1])
            Ψshear[i] = Ψ[i].u[1].coeff*(dΨudy.p(1.0)  + dΨudy.p(-1.0))/2
        end
    end

    println("Making matrices B,A ...")
    y = Polynomial{T, :y}([zero(T), one(T)])
    B   = [innerproduct(Ψ[i], Ψ[j]) for i in 1:m, j in 1:m]
    A1a = [-innerproduct(Ψ[i], y*xderivative(Ψ[j])) for i in 1:m, j in 1:m]
    A1b = [-innerproduct(Ψ[i], vex(Ψ[j])) for i in 1:m, j in 1:m]
    A2  = [innerproduct(Ψ[i], laplacian(Ψ[j])) for i in 1:m, j in 1:m]
    A1 = A1a + A1b

    println("Making quadratic operator N...")
    Ndense = zeros(T, m, m, m)
    for j in 1:m
        print("$j ")
        for k in 1:m
            Ψj_dotgrad_Ψk = dotgrad(Ψ[j], Ψ[k])
            for i in 1:m
                val = -innerproduct(Ψ[i], Ψj_dotgrad_Ψk)
                Ndense[i,j,k] = abs(val) > 1e-15 ? val : zero(T) 
            end
        end
    end
    println()
    N  = SparseBilinear(Ndense)

    # precompute LU decomp of Bf to speed repeated calls to Bf x = b solves
    Bfact = lu(B)
    
    # construct Re-parameterized f and Df functions 
    f(x,R)  = Bfact\(A1*x + (1/R)*(A2*x) + N(x))
    Df(x,R) = Bfact\(A1 + (1/R)*(A2) + derivative(N,x))
    
    ODEModel(α, γ, H, ijkl, Ψ, Ψshear, B, A1, A2, N, Bfact, f, Df)
end

length(model::ODEModel{T}) where T<:Real = length(model.Ψ)
shear(x::Vector{T}, model::ODEModel{T}) where T<:Real = one(T) + dot(x, model.Ψshear)

"""
  `cnab2(x₀, odemodel, R, Δt, Nsteps, nsave)`

Integrate odemodel in time using 2nd-order Crank-Nicolson, Adams-Bashforth time stepping
from initial condition x₀ at Reynolds number R, taking Nsteps steps of length Δt, and
saving every nsave-th time step. Return t and x(t) in vector t, matrix X, where
t[n] = tₙ₋₁ = (n-1)*save*Δt and X[:,n] = x(tₙ₋₁).
"""
function cnab2(x₀, odemodel, R, Δt, Nsteps, nsave; verbosity=0)

    A = model.A1 + (1/R)*model.A2
    N = model.N

    LHS_linear = odemodel.B - (Δt/2)*A
    RHS_linear = odemodel.B + (Δt/2)*A

    LHS_lu = lu(LHS_linear)

    xn  = x₀      # xⁿ
    Nxn1 = N(x₀)  # N(xⁿ⁻¹)
    Nxn  = N(x₀)  # N(xⁿ)

    Nsave = Nsteps÷nsave + 1 # total number of saved steps (+ 1 for step zero)

    X = zeros(length(x₀), Nsave)
    t = (0:Nsave)*Δt

    X[:,1] = x₀

    k = 1
    for n = 1:Nsteps
        verbosity > 0 && n % 10 == 0 && print("$n ")

        Nxn1 .= Nxn    # shift previous N(xⁿ) to N(xⁿ⁻¹)
        Nxn .= N(xn)   # calculate N(xⁿ)
        xn = LHS_lu\(RHS_linear*xn + (3Δt/2)*Nxn - (Δt/2)*Nxn1)

        if n % nsave == 0
            X[:,k] = xn
            k += 1
        end
    end
    t,X
end
