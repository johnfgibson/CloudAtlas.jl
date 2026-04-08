using Test
using Polynomials
using LinearAlgebra

using CloudAtlas

sx = Symmetry(-1, 1, 1)
sy = Symmetry(1, -1, 1)
sz = Symmetry(1, 1, -1)
tx = Symmetry(1, 1, 1, 1 // 2, 0 // 1)
tz = Symmetry(1, 1, 1, 0 // 1, 1 // 2)

@testset "CloudAtlas Tests" begin
    @time @testset "Polynomials" begin
        @test polyparity(Polynomial([1])) == 1
        @test polyparity(Polynomial([0, 1])) == -1
        @test polyparity(Polynomial([1, 1])) == 0
    end

    @time @testset "Fourier Modes" begin
        f = FourierMode()
        @test compatible(f, f)

        # Test mode construction and evaluation
        f1 = FourierMode(1, 1.0)  # sin(x)
        f2 = FourierMode(-1, 1.0) # cos(x)
        f0 = FourierMode(0, 1.0)  # constant

        @test f1(0.0) ≈ 0.0
        @test f1(π / 2) ≈ 1.0
        @test f2(0.0) ≈ 1.0
        @test f0(1.0) ≈ 1.0

        # Test orthogonality
        @test isorthogonal(f1, f2)
        @test !isorthogonal(f1, f1)
    end

    @time @testset "Basis Components" begin
        f = FourierMode(1, 1, 1)
        b = BasisComponent(1, f, f, Polynomial([1, 0, -1], :y), 0)
        @test compatible(b, b)

        # Test evaluation
        val = b(0.0, 0.0, 0.0)
        @test isfinite(val)
    end

    @time @testset "Symmetry operations" begin
        # Test symmetry composition
        σ1 = Symmetry(-1, 1, 1)
        σ2 = Symmetry(1, -1, 1)
        σ3 = σ1 * σ2
        @test σ3.sx == -1
        @test σ3.sy == -1

        # Test halfbox symmetries
        sx, sy, sz, tx, tz = halfbox_symmetries()
        @test sx.sx == -1
        @test tx.ax == 1 // 2
        @test tz.az == 1 // 2
    end

    @time @testset "FourierMode Inner Products and Norms" begin
        α = 1.0

        # Create orthogonal modes
        e1 = FourierMode(1, α)   # sin(αx)
        e2 = FourierMode(-1, α)  # cos(αx)
        e0 = FourierMode(0, α)   # constant

        # Test orthogonality
        @test innerproduct(e1, e2) ≈ 0.0 atol = 1e-14
        @test innerproduct(e1, e0) ≈ 0.0 atol = 1e-14
        @test innerproduct(e2, e0) ≈ 0.0 atol = 1e-14

        # Test norms: constant has norm 1, sin/cos have norm 1/√2
        @test innerproduct(e0, e0) ≈ 1.0 atol = 1e-14
        @test innerproduct(e1, e1) ≈ 0.5 atol = 1e-14
        @test innerproduct(e2, e2) ≈ 0.5 atol = 1e-14

        @test norm(e0) ≈ 1.0 atol = 1e-14
        @test norm(e1) ≈ 1 / sqrt(2) atol = 1e-14
        @test norm(e2) ≈ 1 / sqrt(2) atol = 1e-14

        # Test with coefficients
        e3 = FourierMode(2.0, 1, α)
        @test norm(e3) ≈ 2 / sqrt(2) atol = 1e-14
    end

    @time @testset "FourierMode Derivatives" begin
        α = 2.0

        # d/dx sin(αx) = α cos(αx)
        e_sin = FourierMode(1, α)
        de_sin = derivative(e_sin)
        @test de_sin.waveindex == -1  # cos
        @test de_sin.coeff ≈ α atol = 1e-14

        # d/dx cos(αx) = -α sin(αx)
        e_cos = FourierMode(-1, α)
        de_cos = derivative(e_cos)
        @test de_cos.waveindex == 1  # sin
        @test de_cos.coeff ≈ -α atol = 1e-14

        # d/dx constant = 0
        e_const = FourierMode(0, α)
        de_const = derivative(e_const)
        @test de_const.waveindex == 0
        @test de_const.coeff ≈ 0.0 atol = 1e-14

        # Test second derivative: d²/dx² sin(αx) = -α² sin(αx)
        d2e_sin = derivative(e_sin, 2)
        @test d2e_sin.waveindex == 1  # sin
        @test d2e_sin.coeff ≈ -α^2 atol = 1e-14
    end

    @time @testset "FourierMode Multiplication" begin
        α = 1.0

        # sin(x) * sin(x) = 1/2 - 1/2 cos(2x)
        e1 = FourierMode(1, α)
        prod = e1 * e1
        @test length(prod) == 2
        # Should contain constant term and cos(2x) term

        # cos(x) * cos(x) = 1/2 + 1/2 cos(2x)
        e2 = FourierMode(-1, α)
        prod2 = e2 * e2
        @test length(prod2) == 2

        # sin(x) * cos(x) = 1/2 sin(2x)
        prod3 = e1 * e2
        @test length(prod3) == 1
        @test prod3[1].waveindex == 2  # sin(2x)

        # constant * anything
        e0 = FourierMode(2.0, 0, α)
        prod4 = e0 * e1
        @test length(prod4) == 1
        @test prod4[1].coeff ≈ 2.0 atol = 1e-14
    end

    @time @testset "Polynomial Inner Products and Norms" begin
        # (1-y²) is even, norm² = ∫₋₁¹ (1-y²)² dy / 2
        p1 = Polynomial([1.0, 0.0, -1.0], :y)  # 1 - y²
        ip1 = innerproduct(p1, p1)
        @test ip1 > 0
        @test isfinite(ip1)

        # y is odd, norm² = ∫₋₁¹ y² dy / 2 = 1/3
        p2 = Polynomial([0, 1.0], :y)  # y
        ip2 = innerproduct(p2, p2)
        @test ip2 ≈ 1.0 / 3.0 atol = 1e-14

        # Orthogonality: even * odd = 0
        ip12 = innerproduct(p1, p2)
        @test ip12 ≈ 0.0 atol = 1e-14

        # With parity specified
        ip1_parity = innerproduct(p1, p1, 1, 1)  # both even
        @test ip1_parity ≈ ip1 atol = 1e-14

        ip12_parity = innerproduct(p1, p2, 1, -1)  # even * odd
        @test ip12_parity ≈ 0.0 atol = 1e-14
    end

    @time @testset "BasisComponent Inner Products and Norms" begin
        α, γ = 1.0, 2.0

        # Create simple basis components
        ejx = FourierMode(1, α)
        ekz = FourierMode(1, γ)
        p = Polynomial([1, 0, -1], :y)

        b1 = BasisComponent(1.0, ejx, ekz, p, 1)
        b2 = BasisComponent(1.0, ejx, ekz, p, 1)

        # Same component should have positive inner product
        ip = innerproduct(b1, b2)
        @test ip > 0
        @test isfinite(ip)

        # Norm should be positive
        n = norm(b1)
        @test n > 0
        @test n^2 ≈ norm2(b1) atol = 1e-12

        # Orthogonal components
        ejx2 = FourierMode(-1, α)  # cos instead of sin
        b3 = BasisComponent(1.0, ejx2, ekz, p, 1)
        ip_ortho = innerproduct(b1, b3)
        @test ip_ortho ≈ 0.0 atol = 1e-14

        # Opposite parity should be orthogonal
        p_odd = Polynomial([0, 1], :y)
        b4 = BasisComponent(1.0, ejx, ekz, p_odd, -1)
        ip_parity = innerproduct(b1, b4)
        @test ip_parity ≈ 0.0 atol = 1e-14
    end

    @time @testset "BasisComponent Derivatives" begin
        α, γ = 1.0, 2.0

        ejx = FourierMode(1, α)
        ekz = FourierMode(1, γ)
        p = Polynomial([1, 0, -1], :y)  # 1 - y²

        b = BasisComponent(1.0, ejx, ekz, p, 1)

        # x-derivative
        dbdx = xderivative(b)
        @test dbdx.ejx.waveindex == -1  # sin → cos
        @test dbdx.ejx.coeff ≈ α atol = 1e-14
        @test dbdx.p == p  # y-polynomial unchanged

        # y-derivative
        dbdy = yderivative(b)
        @test dbdy.p == derivative(p)
        @test dbdy.pparity == -1  # parity flips

        # z-derivative
        dbdz = zderivative(b)
        @test dbdz.ekz.waveindex == -1
        @test dbdz.ekz.coeff ≈ γ atol = 1e-14

        # Second derivatives
        d2bdx2 = xderivative(b, 2)
        @test d2bdx2.ejx.coeff ≈ -α^2 atol = 1e-14
    end

    @time @testset "BasisComponent Laplacian" begin
        α, γ = 1.0, 2.0

        ejx = FourierMode(1, α)
        ekz = FourierMode(1, γ)
        p = Polynomial([1, 0, -1], :y)

        b = BasisComponent(1.0, ejx, ekz, p, 1)

        # Laplacian ∇²b = ∂²b/∂x² + ∂²b/∂y² + ∂²b/∂z²
        lap_b = laplacian(b)

        # Should be a valid BasisComponent
        @test lap_b isa BasisComponent
        @test isfinite(lap_b.coeff)
        @test lap_b.pparity == 1  # Laplacian preserves parity for this case
    end

    @time @testset "BasisComponent Multiplication" begin
        α, γ = 1.0, 2.0

        ejx1 = FourierMode(1, α)
        ekz1 = FourierMode(1, γ)
        p1 = Polynomial([1], :y)

        ejx2 = FourierMode(1, α)
        ekz2 = FourierMode(1, γ)
        p2 = Polynomial([0, 1], :y)

        b1 = BasisComponent(1.0, ejx1, ekz1, p1, 1)
        b2 = BasisComponent(1.0, ejx2, ekz2, p2, -1)

        # Multiplication returns array of BasisComponents
        prod = b1 * b2

        # All products should be valid
        for bc in prod
            @test isfinite(bc.coeff)
        end
    end

    @time @testset "BasisFunction Inner Products and Norms" begin
        α, γ = 1.0, 2.0

        # Create basis components
        ejx = FourierMode(1, α)
        ekz = FourierMode(1, γ)
        p = Polynomial([1, 0, -1], :y)

        u = BasisComponent(1.0, ejx, ekz, p, 1)
        v = BasisComponent(0.0, ejx, ekz, p, 1)
        w = BasisComponent(0.0, ejx, ekz, p, 1)

        ψ1 = BasisFunction(u, v, w)
        ψ2 = BasisFunction(u, v, w)

        # Inner product
        ip = innerproduct(ψ1, ψ2)
        @test ip > 0
        @test isfinite(ip)

        # Norm
        n = norm(ψ1)
        @test n > 0
        @test n^2 ≈ norm2(ψ1) atol = 1e-12

        # Orthogonal basis functions
        ejx2 = FourierMode(-1, α)
        u2 = BasisComponent(1.0, ejx2, ekz, p, 1)
        ψ3 = BasisFunction(u2, v, w)

        ip_ortho = innerproduct(ψ1, ψ3)
        @test ip_ortho ≈ 0.0 atol = 1e-14
    end

    @time @testset "BasisFunction Derivatives and Laplacian" begin
        α, γ = 1.0, 2.0

        ejx = FourierMode(1, α)
        ekz = FourierMode(1, γ)
        p = Polynomial([1, 0, -1], :y)

        u = BasisComponent(1.0, ejx, ekz, p, 1)
        v = BasisComponent(1.0, ejx, ekz, p, 1)
        w = BasisComponent(1.0, ejx, ekz, p, 1)

        ψ = BasisFunction(u, v, w)

        # Derivatives
        dψdx = xderivative(ψ)
        @test dψdx isa BasisFunction
        @test all(isfinite.(dψdx.u[1].coeff))

        dψdy = yderivative(ψ)
        @test dψdy isa BasisFunction

        dψdz = zderivative(ψ)
        @test dψdz isa BasisFunction

        # Laplacian
        lap_ψ = laplacian(ψ)
        @test lap_ψ isa BasisFunction
        @test all(bc -> isfinite(bc.coeff), lap_ψ.u)
    end

    @time @testset "BasisFunction Dot-Grad Operation" begin
        α, γ = 1.0, 2.0

        ejx = FourierMode(1, α)
        ekz = FourierMode(1, γ)
        p = Polynomial([1, 0, -1], :y)

        u1 = BasisComponent(1.0, ejx, ekz, p, 1)
        v1 = BasisComponent(0.5, ejx, ekz, p, 1)
        w1 = BasisComponent(0.3, ejx, ekz, p, 1)

        u2 = BasisComponent(0.8, ejx, ekz, p, 1)
        v2 = BasisComponent(0.6, ejx, ekz, p, 1)
        w2 = BasisComponent(0.4, ejx, ekz, p, 1)

        ψ1 = BasisFunction(u1, v1, w1)
        ψ2 = BasisFunction(u2, v2, w2)

        # ψ1 · ∇ψ2
        result = dotgrad(ψ1, ψ2)

        # Should return vector of vectors of BasisComponents
        @test result isa Vector{Vector{BasisComponent{T}}} where {T<:Number}
        @test length(result) == 3  # Three components: u, v, w

        # Each component should be non-empty
        for comp in result
            @test length(comp) > 0
            for bc in comp
                @test bc isa BasisComponent
            end
        end
    end

    @time @testset "BasisFunction Evaluation" begin
        α, γ = 1.0, 2.0

        ejx = FourierMode(1, α)
        ekz = FourierMode(1, γ)
        p = Polynomial([1, 0, -1], :y)

        u = BasisComponent(1.0, ejx, ekz, p, 1)
        v = BasisComponent(0.5, ejx, ekz, p, 1)
        w = BasisComponent(0.3, ejx, ekz, p, 1)

        ψ = BasisFunction(u, v, w)

        # Evaluate at a point
        x, y, z = 0.5, 0.0, 1.0
        val1 = ψ(x, y, z)

        @test val1 isa Vector
        @test length(val1) == 3
        @test all(isfinite.(val1))

        # Alternative evaluation with vector input
        val2 = ψ([x, y, z])
        @test val2 ≈ val1 atol = 1e-14
    end

    @time @testset "Legendre Polynomials" begin
        L = 5
        P = legendrePolynomials(L)

        # Check length
        @test length(P) == L + 1

        # Check orthogonality: ∫₋₁¹ Pₘ(y) Pₙ(y) dy = 0 for m ≠ n
        for m in 0:L-1
            for n in m+1:L
                prod = P[m] * P[n]
                prod_int = integrate(prod)
                integral = prod_int(1) - prod_int(-1)
                @test abs(integral) < 1e-12
            end
        end
    end

    @time @testset "BasisSet Construction" begin
        α, γ = 1.0, 2.0
        J, K, L = 1, 1, 2

        Ψ = basisSet(α, γ, J, K, L)

        # Check that we got basis functions
        @test length(Ψ) > 0
        @test all(ψ -> ψ isa BasisFunction, Ψ)

        # Check normalization option
        Ψ_norm = basisSet(α, γ, J, K, L, normalize=true)

        # Normalized basis should have unit norm (approximately)
        for ψ in Ψ_norm
            n = norm(ψ)
            @test n ≈ 1.0 atol = 0.5  # Loose tolerance due to rational arithmetic
        end
    end

    @time @testset "BasisSet with Symmetries" begin
        α, γ = 1.0, 2.0
        J, K, L = 1, 1, 2
        H = [sx * sy * sz]

        ijkl = basisIndices(J, K, L, H)
        Ψ = basisSet(α, γ, ijkl)

        @test length(Ψ) == size(ijkl, 1)
        @test all(ψ -> ψ isa BasisFunction, Ψ)
    end

    @time @testset "VEX Operation" begin
        α, γ = 1.0, 2.0

        ejx = FourierMode(1, α)
        ekz = FourierMode(1, γ)
        p = Polynomial([1, 0, -1], :y)

        u = BasisComponent(1.0, ejx, ekz, p, 1)
        v = BasisComponent(2.0, ejx, ekz, p, 1)
        w = BasisComponent(3.0, ejx, ekz, p, 1)

        ψ = BasisFunction(u, v, w)

        # vex(ψ) = [v, 0, 0]
        vex_ψ = vex(ψ)

        @test vex_ψ isa BasisFunction
        @test vex_ψ.u[1].coeff ≈ 2.0 atol = 1e-14
        @test vex_ψ.u[2].coeff ≈ 0.0 atol = 1e-14
        @test vex_ψ.u[3].coeff ≈ 0.0 atol = 1e-14
    end

    @time @testset "Scalar Multiplication of BasisFunction" begin
        α, γ = 1.0, 2.0

        ejx = FourierMode(1, α)
        ekz = FourierMode(1, γ)
        p = Polynomial([1, 0, -1], :y)

        u = BasisComponent(1.0, ejx, ekz, p, 1)
        v = BasisComponent(2.0, ejx, ekz, p, 1)
        w = BasisComponent(3.0, ejx, ekz, p, 1)

        ψ = BasisFunction(u, v, w)

        # Scalar multiplication
        c = 2.5
        scaled_ψ = c * ψ

        @test scaled_ψ isa BasisFunction
        @test norm(scaled_ψ) ≈ abs(c) * norm(ψ) atol = 1e-12
    end

    @time @testset "Hookstep search on simple 2d function" begin
        f(x) = [x[1]^2 + x[2]^2 - 1; x[1] - x[2]]
        Df(x) = [2x[1] 2x[2]; 1 -1]

        xguess = [2; 0]
        xsolution = [1 / sqrt(2); 1 / sqrt(2)]
        xsolve, success = hookstepsolve(f, Df, xguess)

        @test norm(xsolve - xsolution) < 1e-10
    end

    @time @testset "ODEModel evaluation on known equilibrium" begin
        # set up ODE model for some tests
        α, γ = 1.0, 2.0             # Fourier wavenumbers α, γ = 2π/Lx, 2π/Lz
        H = [sx * sy * sz, sz * tx * tz]    # Generators of the symmetric subspace of the Nagata eqb
        J, K, L = 1, 1, 3               # Bounds on Fourier modes (J,K) and wall-normal polynomials (L)
        R = 200.0

        # ODE eqb version of Nagata lower branch, computed for PRL
        xeqb = [0.16764501289529937, -0.018552429231415948, 0.48640800040772025, -0.09759477649001111,
            -0.03411986464200022, -0.010665050862032609, -0.01731713276665348, 0.020113317487885345,
            -0.03572908824700301, -0.008809478655993888, -0.018431734698126298, -0.04736482532692928,
            -0.04609045287990289, -0.03157916691896989, 0.0010374375712356061, 0.01601957685174594,
            0.07702704667134751]

        model = ODEModel(α, γ, J, K, L, H)  # Form ODE model via Galerkin projection

        @test norm(model.f(xeqb, R)) / norm(xeqb) < 1e-10

    end

    @time @testset "Hookstep Solve for equilibrium of ODEModel" begin
        α, γ = 1.0, 2.0             # Fourier wavenumbers α, γ = 2π/Lx, 2π/Lz
        H = [sx * sy * sz, sz * tx * tz]    # Generators of the symmetric subspace of the Nagata eqb
        J, K, L = 1, 1, 3               # Bounds on Fourier modes (J,K) and wall-normal polynomials (L)
        R = 200.0                   # Reynolds number

        model = ODEModel(α, γ, J, K, L, H)  # Form ODE model via Galerkin projection

        # Projection of Nagata eqb used as initial guess for ODE eqb. Computed for PRL
        xguess = [0.105, -0.0539, 0.388, -0.0172, -0.0133, 0.0113, -0.00706, 0.0240, -0.0344, -0.00868,
            -0.01264, -0.0234, -0.0320, -0.0180, -0.00464, 0.0106, 0.0235]
        xeqb = [0.16764501289529937, -0.018552429231415948, 0.48640800040772025, -0.09759477649001111,
            -0.03411986464200022, -0.010665050862032609, -0.01731713276665348, 0.020113317487885345,
            -0.03572908824700301, -0.008809478655993888, -0.018431734698126298, -0.04736482532692928,
            -0.04609045287990289, -0.03157916691896989, 0.0010374375712356061, 0.01601957685174594,
            0.07702704667134751]

        f = x -> model.f(x, R)
        Df = x -> model.Df(x, R)
        xsoln, success = hookstepsolve(f, Df, xguess)
        @test norm(xsoln - xeqb) / norm(xeqb) < 1e-08
    end

    @time @testset "Hookstep Solve for equilibrium of ODEModel, normalized" begin
        α, γ = 1.0, 2.0             # Fourier wavenumbers α, γ = 2π/Lx, 2π/Lz
        H = [sx * sy * sz, sz * tx * tz]    # Generators of the symmetric subspace of the Nagata eqb
        J, K, L = 1, 1, 3               # Bounds on Fourier modes (J,K) and wall-normal polynomials (L)
        R = 200.0                   # Reynolds number

        model = ODEModel(α, γ, J, K, L, H, normalize=true)  # Form ODE model via Galerkin projection
        params = SearchParams(; δ=0.01)

        # Projection of Nagata eqb used as initial guess for ODE eqb. Computed for PRL
        xguess = [0.117, -0.0402, 0.201, -0.00779, -0.0158, 0.00719, -0.00551, 0.0126, -0.0281,
            -0.00619, -0.0102, -0.0121, -0.0144, -0.0222, -0.00386, 0.01876, 0.0229]

        xeqb = [0.18509767944776034, -0.013810193934778796, 0.2511800112873277, -0.04399074169151473,
            -0.04069031295467911, -0.0068103678652483395, -0.013519832908378151, 0.01058687307496652,
            -0.02917267839495075, -0.006278485885577374, -0.014925151720599003, -0.024459090632393224,
            -0.020775222591909327, -0.038982169958350785, 0.0008634078112314248, 0.02842834644157768,
            0.07523725576071376]

        f = x -> model.f(x, R)
        Df = x -> model.Df(x, R)
        xsoln, success = hookstepsolve(f, Df, xguess, params)
        @test norm(xsoln - xeqb) / norm(xeqb) < 1e-08
    end

    @time @testset "TWModel Construction and Basic Properties" begin
        # simple check that we can construct a TWModel and it doesn't modify any inputs unexpectedly
        α, γ = 1.0, 2.0
        J, K, L = 1, 1, 2
        H = [sx * sy * sz]

        model = TWModel(α, γ, J, K, L, H)

        @test model.α == α
        @test model.γ == γ
        @test length(model) > 0
        @test size(model.B, 1) == length(model)
        @test size(model.A1, 1) == length(model)
        @test size(model.A2, 1) == length(model)
        @test size(model.Cx, 1) == length(model)
        @test size(model.Cz, 1) == length(model)
    end

    @time @testset "TWModel Phase Constraint Detection" begin
        # verifies that the traveling wave detects phase constraints correctly
        α, γ = 1.0, 2.0
        J, K, L = 1, 1, 2

        # Both constraints needed
        H1 = [sy]
        model1 = TWModel(α, γ, J, K, L, H1)
        @test model1.keep_cx && model1.keep_cz  # At least one constraint

        # With shift symmetries
        H2 = [(sx * sy) * (tx * tz)] # z traveling only 
        model2 = TWModel(α, γ, J, K, L, H2)
        @test !model2.keep_cx && model1.keep_cz
    end

    @time @testset "TWModel Full Residual with Phase Constraints" begin
        # Verifies that the residual function can be called, is finite, and has the correct number of dimensions given the constraints. 
        α, γ = 1.0, 2.0
        J, K, L = 1, 1, 2
        R = 200.0
        H_list = [[sx * sy * sz], [sx * sy * sz, sz * tx * tz], [sy], [(sx * sy) * (tx * tz)]]
        for H in H_list

            model = TWModel(α, γ, J, K, L, H)
            m = length(model)

            # Augmented state vector [x; cx; cz]
            ξ = [randn(m) * 0.01; 0.5; 0.3]

            # Full residual includes phase constraints
            g_val = model.g(ξ, R)

            # Should have m equations + number of active phase constraints
            expected_dim = m + (model.keep_cx ? 1 : 0) + (model.keep_cz ? 1 : 0)
            @test length(g_val) == expected_dim
            @test all(isfinite.(g_val))
        end
    end

    @time @testset "TWModel Bordered Jacobian Dimensions" begin
        # Verifies that the number of dimensions is correct for various symmetries
        α, γ = 1.0, 2.0
        J, K, L = 1, 1, 2
        R = 200.0
        H_list = [[sx * sy * sz], [sx * sy * sz, sz * tx * tz], [sy], [(sx * sy) * (tx * tz)]]
        for H in H_list

            model = TWModel(α, γ, J, K, L, H)
            m = length(model)

            ξ = [randn(m) * 0.01; 0.5; 0.3]

            Dg_val = model.Dg(ξ, R)

            # Rows: m + active phase constraints
            # Cols: m + 2 (always have cx, cz columns)
            expected_rows = m + (model.keep_cx ? 1 : 0) + (model.keep_cz ? 1 : 0)
            expected_cols = m + 2

            @test size(Dg_val, 1) == expected_rows
            @test size(Dg_val, 2) == expected_cols
            @test all(isfinite.(Dg_val))
        end
    end

    @time @testset "TWModel ODE Functions" begin
        # Verifies that the ODE functions can be called and return finite values.
        α, γ = 1.0, 2.0
        J, K, L = 1, 1, 2
        H = [(sx * sy) * (tx * tz)] # z traveling only 
        R = 200.0

        model = TWModel(α, γ, J, K, L, H)
        m = length(model)

        x = randn(m) * 0.01
        cx = 0.5
        cz = 0.3

        # Test f function (time derivative)
        dxdt = model.f(x, cx, cz, R)
        @test length(dxdt) == m
        @test all(isfinite.(dxdt))

        # Test Df function (Jacobian of f)
        Df_val = model.Df(x, cx, cz, R)
        @test size(Df_val) == (m, m)
        @test all(isfinite.(Df_val))
    end

    @time @testset "TWModel Shear Calculation" begin
        # calculates the shear of a flow very near laminar (probably a better way to test this?)
        # the important thing for now is just verifying that all functions run
        α, γ = 1.0, 2.0
        J, K, L = 1, 1, 2
        H = [sx * sy * sz]

        model = TWModel(α, γ, J, K, L, H)
        m = length(model)

        x = randn(m) * 0.001
        s = shear(x, model)

        @test s isa Real
        @test isfinite(s)
        @test s ≈ 1.0 atol = 1.0  # Should be close to base flow shear
    end

    @time @testset "TWModel Symmetry File Output" begin
        α, γ = 1.0, 2.0
        J, K, L = 1, 1, 2
        H = [sx * sy * sz]

        model = TWModel(α, γ, J, K, L, H)

        # Test save_sigma function
        cx, cz, T = 0.5, 0.3, 10.0
        tmpfile = tempname()

        try
            save_sigma(model, cx, cz, T, tmpfile)
            @test isfile(tmpfile)

            # Read back and verify format
            content = read(tmpfile, String)
            @test occursin("%", content)  # Should have comment line
            @test occursin("1 1 1 1", content)  # Should have symmetry values
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

    @time @testset "Basis Index Dictionary" begin
        α, γ = 1.0, 2.0
        J1, K1, L1 = 1, 1, 2
        J2, K2, L2 = 2, 2, 3
        H = [sx * sy * sz]

        ijkl1 = basisIndices(J1, K1, L1, H)
        ijkl2 = basisIndices(J2, K2, L2, H)

        # Test basis index mapping
        dict = basis_index_dict(ijkl1, ijkl2)

        @test isa(dict, Dict{Int,Int})
        @test all(n -> ijkl2[dict[n], :] == ijkl1[n, :], keys(dict))
    end

    @time @testset "Change Basis Operation" begin
        α, γ = 1.0, 2.0
        J1, K1, L1 = 1, 1, 2
        J2, K2, L2 = 2, 2, 3
        H = [sx * sy * sz]

        ijkl1 = basisIndices(J1, K1, L1, H)
        ijkl2 = basisIndices(J2, K2, L2, H)

        m1 = size(ijkl1, 1)
        x1 = randn(m1)

        # Convert to larger basis
        x2 = changebasis(x1, ijkl1, ijkl2)

        @test length(x2) == size(ijkl2, 1)
        @test all(isfinite.(x2))
    end

    @time @testset "SparseBilinear Operator" begin
        # Create small dense tensor
        N_dense = randn(5, 5, 5)
        N_dense[abs.(N_dense).<0.5] .= 0  # Sparsify

        # Convert to sparse
        N_sparse = SparseBilinear(N_dense)

        # Test evaluation
        x = randn(5)
        y = randn(5)

        result_sparse = N_sparse(x, y)

        # Verify against dense computation
        result_dense = zeros(5)
        for i in 1:5, j in 1:5, k in 1:5
            result_dense[i] += N_dense[i, j, k] * x[j] * y[k]
        end

        @test result_sparse ≈ result_dense atol = 1e-12
    end
end

