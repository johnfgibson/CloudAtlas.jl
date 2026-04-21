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
    #Ψshear::Vector{T}            # value of dΨudy at walls, used to calculate shear rate
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
    
    ODEModel(α, γ, H, ijkl, Ψ, B, A1, A2, N, Bfact, f, Df)
end

length(model::ODEModel{T}) where T<:Real = length(model.Ψ)

function shear_function(Ψ::Vector{BasisFunction{T}}) where T<:Real
    m = length(Ψ)
    Ψshear = zeros(T, m) # vector of shear rates of each Ψ[j]
    for i=1:m
        if Ψ[i].u[1].ejx.waveindex == 0 && Ψ[i].u[1].ekz.waveindex == 0
            dΨudy = yderivative(Ψ[i].u[1])
            Ψshear[i] = Ψ[i].u[1].coeff*(dΨudy.p(1.0)  + dΨudy.p(-1.0))/2
        end
    end

    # build and return function shear(x) = 1 + sum x[j] shear(Ψ[j])
    shear(x::Vector{T}) = one(T) + dot(x, Ψshear)
    return shear 
end


"""
    dissipation_term(f, g)

Compute D_ij = < curl(f), curl(g) > for two basis functions.
Expands the curl components to avoid constructing summed BasisComponents.
"""
function dissipation_term(f::BasisFunction{T}, g::BasisFunction{T}) where {T<:Real}
    # f.u[1] = u, f.u[2] = v, f.u[3] = w
    f_dy_u = yderivative(f.u[1])
    f_dz_u = zderivative(f.u[1])
    f_dx_v = xderivative(f.u[2])
    f_dz_v = zderivative(f.u[2])
    f_dx_w = xderivative(f.u[3])
    f_dy_w = yderivative(f.u[3])

    g_dy_u = yderivative(g.u[1])
    g_dz_u = zderivative(g.u[1])
    g_dx_v = xderivative(g.u[2])
    g_dz_v = zderivative(g.u[2])
    g_dx_w = xderivative(g.u[3])
    g_dy_w = yderivative(g.u[3])

    # X-component: <dy w - dz v, dy w - dz v>
    val_x = innerproduct(f_dy_w, g_dy_w) -
            innerproduct(f_dy_w, g_dz_v) -
            innerproduct(f_dz_v, g_dy_w) +
            innerproduct(f_dz_v, g_dz_v)

    # Y-component: <dz u - dx w, dz u - dx w>
    val_y = innerproduct(f_dz_u, g_dz_u) -
            innerproduct(f_dz_u, g_dx_w) -
            innerproduct(f_dx_w, g_dz_u) +
            innerproduct(f_dx_w, g_dx_w)

    # Z-component: <dx v - dy u, dx v - dy u>
    val_z = innerproduct(f_dx_v, g_dx_v) -
            innerproduct(f_dx_v, g_dy_u) -
            innerproduct(f_dy_u, g_dx_v) +
            innerproduct(f_dy_u, g_dy_u)

    return val_x + val_y + val_z
end

function dissipation_function(Ψ::Vector{BasisFunction{T}}) where T<:Real
    m = length(Ψ)
    D = zeros(T, m, m)
    for j in 1:m
        for i in 1:j
            val = dissipation_term(Ψ[i], Ψ[j])
            D[i, j] = val
            D[j, i] = val
        end
    end
    # build and return function dissipation(x) = 1 + x D x
    dissipation(x::Vector{T}) = one(T) + dot(x, D*x)
    return dissipation
end

"""
  `cnab2(x₀, odemodel, R, Δt, Nsteps, nsave)`

Integrate odemodel in time using 2nd-order Crank-Nicolson, Adams-Bashforth time stepping
from initial condition x₀ at Reynolds number R, taking Nsteps steps of length Δt, and
saving every nsave-th time step. Return t and x(t) in vector t, matrix X, where
t[n] = tₙ₋₁ = (n-1)*save*Δt and X[:,n] = x(tₙ₋₁).
"""
function cnab2(x₀, odemodel, R, Δt, Nsteps, nsave; verbosity=0)

    A = odemodel.A1 + (1/R)*odemodel.A2
    N = odemodel.N

    LHS_linear = odemodel.B - (Δt/2)*A
    RHS_linear = odemodel.B + (Δt/2)*A

    LHS_lu = lu(LHS_linear)

    xn  = x₀      # xⁿ
    Nxn1 = N(x₀)  # N(xⁿ⁻¹)
    Nxn  = N(x₀)  # N(xⁿ)

    Nsave = Nsteps÷nsave # total number of saved steps (+ 1 for step zero)

    X = zeros(Nsave, length(x₀))
    t = (0:Nsave-1)*Δt

    X[1,:] = x₀

    k = 1
    for n = 1:Nsteps
        verbosity > 0 && n % 10 == 0 && print("$n ")
        Nxn1 .= Nxn    # shift previous N(xⁿ) to N(xⁿ⁻¹)
        Nxn .= N(xn)   # calculate N(xⁿ)
        xn = LHS_lu\(RHS_linear*xn + (3Δt/2)*Nxn - (Δt/2)*Nxn1)

        if n % nsave == 0
            X[k,:] = xn
            k += 1
        end
    end
    t,X
end
