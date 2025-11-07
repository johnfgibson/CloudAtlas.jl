using Test

using CloudAtlas
using Polynomials
using LinearAlgebra

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

    @time @testset "Hookstep Solve" begin
        α, γ = 1.0, 2.0                     # Fourier wavenumbers α, γ = 2π/Lx, 2π/Lz
        J,K,L = 1,1,3                       # Bounds on Fourier modes (J,K) and wall-normal polynomials (L)
        R = 250
        guessNorm = 0.2

        id = Symmetry()
        sx = Symmetry(-1,1,1)
        sy = Symmetry(1,-1,1)
        sz = Symmetry(1,1,-1)
        tx = Symmetry(1,1,1, 1//2, 0//1)
        tz = Symmetry(1,1,1, 0//1, 1//2)

        H = [sx*sy*sz, sz*tx*tz]            # Generators of the symmetric subspace of the Nagata eqb

        ijkl = basisIndices(J,K,L, H) # Compute index set of H-symmetric basis elements Ψijkl
        Ψ = basisSet(α, γ, ijkl)      # Compute basis elements Ψijkl in the index set
        f, Df = ODEModel(Ψ)         # Do Galerkin projection, return f(x,R) for ODE dx/dt = f(x,R)

        Nmodes = length(Ψ) 
        xguess = normalize(randn(Nmodes)) * guessNorm

        xsoln, success = hookstepsolve(f, Df, xguess, ftol=1e-08, xtol=1e-12, Nnewton=30,Nhook=8,verbosity=0)

        @test success
    end
end
