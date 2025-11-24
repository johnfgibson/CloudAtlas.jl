"""
Traveling Wave Model for plane Couette flow
"""
struct TWModel{T<:Real,TB,TF,TDF,TG,TDG}
    α::T                         # streamwise wavenumber α = 2pi/Lx
    γ::T                         # spanwise wavenumber γ = 2pi/Lz
    H::Vector{Symmetry}          # symmetry subgroup
    ijkl::Matrix{Int}            # set of allowed i,j,k,l indices 
    Ψ::Vector{BasisFunction{T}}  # basis set built from allowed i,j,k,l
    Ψshear::Vector{T}            # value of dΨudy at walls, used to calculate shear rate
    B::AbstractMatrix{T}         # inner product matrix 
    A1::AbstractMatrix{T}        # linear term arising from base flow
    A2::AbstractMatrix{T}        # linear term arising from Laplacian
    Cx::AbstractMatrix{T}        # x-translation operator (linear term from cx*∂/∂x)
    Cz::AbstractMatrix{T}        # z-translation operator (linear term from cz*∂/∂z)
    N::SparseBilinear{T}         # nonlinear term from u dot grad u
    Bfact::TB                    # LU factorization of inner product matrix
    m::Int                       # dimension of state space (# of basis functions)
    keep_cx::Bool                # whether to enforce x-phase constraint
    keep_cz::Bool                # whether to enforce z-phase constraint
    f::TF                        # RHS function for ODE: dx/dt = f(x, cx, cz, R)
    Df::TDF                      # Df(x, cx, cz, R), Jacobian w.r.t. x only
    g::TG                        # Residual function: g(ξ, R) = 0 where ξ = [x; cx; cz]
    Dg::TDG                      # Dg(ξ, R), full bordered Jacobian
end

"""
    TWModel(α, γ, J, K, L, H; normalize=false)

Construct a Traveling Wave model for plane Couette flow.

The model solves for traveling wave solutions to:
    c_x ∂u/∂x + c_z ∂u/∂z + u·∇u + v e_x + ∇p - ∇²u = 0

where (cx, cz) is the wave speed and the solution satisfies:
    u(x+cx*T, y, z+cz*T, t+T) = u(x, y, z, t)

Phase constraints are automatically determined from the symmetry group H:
- If H contains a non-trivial x-shift, no x-phase constraint is needed
- If H contains a non-trivial z-shift, no z-phase constraint is needed
"""
function TWModel(
    α::T,
    γ::T,
    J::Int,
    K::Int,
    L::Int,
    H::Vector{Symmetry};
    normalize=false,
    tol_phase=1e-12,
) where {T<:Real}
    ijkl = basisIndices(J, K, L, H)
    m = size(ijkl, 1)
    println("J,K,L,m == $J,$K,$L,$m")
    println("(2J+1)(2K+1)(2L+1) + 1 == $((2J+1)*(2K+1)*(2L+1) + 1)")

    Ψ = basisSet(α, γ, ijkl; normalize=normalize)

    # Compute shear values at walls
    Ψshear = zeros(T, m)
    for i in 1:m
        if Ψ[i].u[1].ejx.waveindex == 0 && Ψ[i].u[1].ekz.waveindex == 0
            dΨudy = yderivative(Ψ[i].u[1])
            Ψshear[i] = Ψ[i].u[1].coeff * (dΨudy.p(1.0) + dΨudy.p(-1.0)) / 2
        end
    end

    println("Making matrices B, A1, A2, Cx, Cz...")
    y = Polynomial{T,:y}([zero(T), one(T)])
    B = [innerproduct(Ψ[i], Ψ[j]) for i in 1:m, j in 1:m]
    A1a = [-innerproduct(Ψ[i], y * xderivative(Ψ[j])) for i in 1:m, j in 1:m]
    A1b = [-innerproduct(Ψ[i], vex(Ψ[j])) for i in 1:m, j in 1:m]
    A2 = [innerproduct(Ψ[i], laplacian(Ψ[j])) for i in 1:m, j in 1:m]
    Cx = [innerproduct(Ψ[i], xderivative(Ψ[j])) for i in 1:m, j in 1:m]
    Cz = [innerproduct(Ψ[i], zderivative(Ψ[j])) for i in 1:m, j in 1:m]
    A1 = A1a + A1b

    # Determine which constraints to keep
    # There should be a better way to do this based on H, but 
    # decided that this has less edgecases and is therefore cleaner 
    cx_nonzero = any(abs.(Cx) .> tol_phase)
    cz_nonzero = any(abs.(Cz) .> tol_phase)

    # Need phase constraint only if the corresponding operator is non-zero
    keep_cx = cx_nonzero
    keep_cz = cz_nonzero

    println("Phase constraints: keep_cx = $keep_cx, keep_cz = $keep_cz")

    println("Making quadratic operator N...")
    Ndense = zeros(T, m, m, m)
    for j in 1:m
        print("$j ")
        for k in 1:m
            Ψj_dotgrad_Ψk = dotgrad(Ψ[j], Ψ[k])
            for i in 1:m
                val = -innerproduct(Ψ[i], Ψj_dotgrad_Ψk)
                Ndense[i, j, k] = abs(val) > 1e-15 ? val : zero(T)
            end
        end
    end
    println()
    N = SparseBilinear(Ndense)

    Bfact = lu(B)

    # ODE form: dx/dt = f(x, cx, cz, R)
    # This is the "lab frame" equation: ∂u/∂t = -u·∇u - v ex - cx ∂u/∂x - cz ∂u/∂z - ∇p + (1/R)∇²u
    function f(x::AbstractVector, cx::Real, cz::Real, R::Real)
        return Bfact \ (A1 * x + (1 / R) * (A2 * x) + cx * Cx * x + cz * Cz * x + N(x))
    end

    # Jacobian of f with respect to x only
    function Df(x::AbstractVector, cx::Real, cz::Real, R::Real)
        return Bfact \ (A1 + (1 / R) * A2 + cx * Cx + cz * Cz + derivative(N, x))
    end

    # Residual function for traveling wave: g(ξ, R) = 0
    # where ξ = [x; cx; cz] is the augmented state vector
    function g(ξ::AbstractVector, R::Real)
        x = ξ[1:m]
        cx = ξ[m + 1]
        cz = ξ[m + 2]

        # Main residual: (cx*Cx + cz*Cz)*x - A1*x - (1/R)*A2*x - N(x) = 0
        # This comes from setting dx/dt = 0 in the moving frame
        residual = (cx * Cx + cz * Cz) * x - A1 * x - (1 / R) * A2 * x - N(x)

        # Dimension of full system
        dim = m + (keep_cx ? 1 : 0) + (keep_cz ? 1 : 0)
        g_full = zeros(eltype(residual), dim)
        g_full[1:m] = residual

        # Phase constraints (if needed)
        idx = m
        if keep_cx
            idx += 1
            g_full[idx] = dot(Cx * x, x)  # x' * Cx * x = 0 fixes x-phase
        end
        if keep_cz
            idx += 1
            g_full[idx] = dot(Cz * x, x)  # x' * Cz * x = 0 fixes z-phase
        end

        return g_full
    end

    # Bordered Jacobian: Dg(ξ, R)
    function Dg(ξ::AbstractVector, R::Real)
        x = ξ[1:m]
        cx = ξ[m + 1]
        cz = ξ[m + 2]

        # Jacobian of main residual with respect to x
        Jx = cx * Cx + cz * Cz - A1 - (1 / R) * A2 - derivative(N, x)

        # Partial derivatives with respect to cx and cz
        ∂r_∂cx = Cx * x
        ∂r_∂cz = Cz * x

        # Build bordered Jacobian
        m_rows = m + (keep_cx ? 1 : 0) + (keep_cz ? 1 : 0)
        m_cols = m + 2  # Always have cx, cz columns

        Dg_matrix = zeros(eltype(Jx), m_rows, m_cols)
        Dg_matrix[1:m, 1:m] = Jx
        Dg_matrix[1:m, m + 1] = ∂r_∂cx
        Dg_matrix[1:m, m + 2] = ∂r_∂cz

        # Add phase constraint rows
        row = m
        if keep_cx
            row += 1
            Dg_matrix[row, 1:m] = 2 * (Cx * x)'  # d/dx(x'*Cx*x) = 2*x'*Cx
            # ∂/∂cx and ∂/∂cz of the phase constraint are zero
        end
        if keep_cz
            row += 1
            Dg_matrix[row, 1:m] = 2 * (Cz * x)'
        end

        return Dg_matrix
    end

    return TWModel(
        α,
        γ,
        H,
        ijkl,
        Ψ,
        Ψshear,
        B,
        A1,
        A2,
        Cx,
        Cz,
        N,
        Bfact,
        m,
        keep_cx,
        keep_cz,
        f,
        Df,
        g,
        Dg,
    )
end

"""
    residual(model::TWModel, x, cx, cz, R)

Compute the residual of the traveling wave equation without phase constraints.
Returns just the m-dimensional residual vector.
"""
function residual(model::TWModel, x::AbstractVector, cx::Real, cz::Real, R::Real)
    return (cx * model.Cx + cz * model.Cz) * x - model.A1 * x - (1 / R) * model.A2 * x -
           model.N(x)
end

"""
    jacobian(model::TWModel, x, cx, cz, R)

Compute the Jacobian of the residual with respect to x only.
"""
function jacobian(model::TWModel, x::AbstractVector, cx::Real, cz::Real, R::Real)
    return cx * model.Cx + cz * model.Cz - model.A1 - (1 / R) * model.A2 -
           derivative(model.N, x)
end

"""
    save_sigma(model::TWModel, cx, cz, T, filename)

Save symmetry file in Channelflow format.
"""
function save_sigma(model::TWModel, cx::Real, cz::Real, T::Real, filename::String)
    open(filename, "w") do file
        az = cz ≈ 0 ? 0 : round(-cz * T; sigdigits=6)
        ax = cx ≈ 0 ? 0 : round(-cx * T; sigdigits=6)
        write(file, "% 1\n1 1 1 1 $(ax) $(az)")
    end
end

length(model::TWModel{T}) where {T<:Real} = model.m
shear(x::Vector{T}, model::TWModel{T}) where {T<:Real} = one(T) + dot(x, model.Ψshear)
