using Test

using CloudAtlas
using Polynomials
using LinearAlgebra

sx = Symmetry(-1, 1, 1)
sy = Symmetry(1, -1, 1)
sz = Symmetry(1, 1, -1)
tx = Symmetry(1, 1, 1, 1//2, 0//1)
tz = Symmetry(1, 1, 1, 0//1, 1//2)

@testset "CloudAtlas Tests" begin
    @time @testset "Polynomials" begin
        @test polyparity(Polynomial([1])) == 1
    end

    @time @testset "Fourier Modes" begin
        f = FourierMode()
        @test compatible(f, f)
    end

    @time @testset "Basis Components" begin
        f = FourierMode(1, 1, 1)
        b = BasisComponent(1, f, f, Polynomial([1, 0, -1], :y), 0)
        @test compatible(b, b)
    end

    @time @testset "ODEModel evaluation on known equilibrium" begin

        # set up ODE model for some tests
        α, γ = 1.0, 2.0                     # Fourier wavenumbers α, γ = 2π/Lx, 2π/Lz
        J,K,L = 1,1,3                       # Bounds on Fourier modes (J,K) and wall-normal polynomials (L)
        R = 200.0
        
        # ODE eqb version of Nagata lower branch, computed for PRL
        xeqb = [0.16764501289529937, -0.018552429231415948, 0.48640800040772025, -0.09759477649001111,
                -0.03411986464200022, -0.010665050862032609, -0.01731713276665348, 0.020113317487885345,
                -0.03572908824700301, -0.008809478655993888, -0.018431734698126298, -0.04736482532692928,
                -0.04609045287990289, -0.03157916691896989, 0.0010374375712356061, 0.01601957685174594,
                 0.07702704667134751]
        
        H = [sx * sy * sz, sz * tx * tz]  # Generators of the symmetric subspace of the Nagata eqb

        ijkl = basisIndices(J, K, L, H)   # Compute index set of H-symmetric basis elements Ψijkl
        Ψ = basisSet(α, γ, ijkl)          # Compute basis elements Ψijkl in the index set
        model = ODEModel(Ψ)               # Do Galerkin projection, return f(x,R) for ODE dx/dt = f(x,R)
        
        @test norm(model.f(xeqb, R))/norm(xeqb) < 1e-10
        
    end

    @time @testset "Hookstep Solve" begin
        α, γ = 1.0, 2.0                     # Fourier wavenumbers α, γ = 2π/Lx, 2π/Lz
        J,K,L = 1,1,3                       # Bounds on Fourier modes (J,K) and wall-normal polynomials (L)
        guessNorm = 0.2
        R = 200.0

        H = [sx * sy * sz, sz * tx * tz]    # Generators of the symmetric subspace of the Nagata eqb

        ijkl = basisIndices(J, K, L, H)     # Compute index set of H-symmetric basis elements Ψijkl
        Ψ = basisSet(α, γ, ijkl)            # Compute basis elements Ψijkl in the index set
        model = ODEModel(Ψ)                 # Do Galerkin projection, return f(x,R) for ODE dx/dt = f(x,R)
        Nmodes = length(Ψ)

        # Projection of Nagata eqb used as initial guess for ODE eqb. Computed for PRL
        xguess = [0.105, -0.0539, 0.388, -0.0172, -0.0133, 0.0113, -0.00706, 0.0240, -0.0344, -0.00868,
                  -0.01264, -0.0234, -0.0320, -0.0180, -0.00464, 0.0106, 0.0235]
        xeqb = [0.16764501289529937, -0.018552429231415948, 0.48640800040772025, -0.09759477649001111,
                -0.03411986464200022, -0.010665050862032609, -0.01731713276665348, 0.020113317487885345,
                -0.03572908824700301, -0.008809478655993888, -0.018431734698126298, -0.04736482532692928,
                -0.04609045287990289, -0.03157916691896989, 0.0010374375712356061, 0.01601957685174594,
                 0.07702704667134751]

        fr = x -> model.f(x,R)
        Dfr = x -> model.Df(x,R)
        xsoln, success = hookstepsolve(fr, Dfr, xguess)
        @test norm(xsoln - xeqb)/norm(xeqb) < 1e-08
    end

    @time @testset "Hookstep Solve with ODEModel and SearchParams" begin
        α, γ = 1.0, 2.0                     # Fourier wavenumbers α, γ = 2π/Lx, 2π/Lz
        J,K,L = 1,1,3                       # Bounds on Fourier modes (J,K) and wall-normal polynomials (L)

        params = SearchParams(; R = 200.0)

        H = [sx * sy * sz, sz * tx * tz]    # Generators of the symmetric subspace of the Nagata eqb

        ijkl = basisIndices(J, K, L, H)     # Compute index set of H-symmetric basis elements Ψijkl
        Ψ = basisSet(α, γ, ijkl)            # Compute basis elements Ψijkl in the index set
        model = ODEModel(Ψ)                 # Do Galerkin projection, return f(x,R) for ODE dx/dt = f(x,R)
        Nmodes = length(Ψ)

        # Projection of Nagata eqb used as initial guess for ODE eqb. Computed for PRL
        xguess = [0.105, -0.0539, 0.388, -0.0172, -0.0133, 0.0113, -0.00706, 0.0240, -0.0344, -0.00868,
                  -0.01264, -0.0234, -0.0320, -0.0180, -0.00464, 0.0106, 0.0235]
        xeqb = [0.16764501289529937, -0.018552429231415948, 0.48640800040772025, -0.09759477649001111,
                -0.03411986464200022, -0.010665050862032609, -0.01731713276665348, 0.020113317487885345,
                -0.03572908824700301, -0.008809478655993888, -0.018431734698126298, -0.04736482532692928,
                -0.04609045287990289, -0.03157916691896989, 0.0010374375712356061, 0.01601957685174594,
                 0.07702704667134751]

        xsoln, success = hookstepsolve(model, xguess, params)
        @test norm(xsoln - xeqb)/norm(xeqb) < 1e-08
    end
end
