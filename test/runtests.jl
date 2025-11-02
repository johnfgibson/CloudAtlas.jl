using Test

using CloudAtlas
using Polynomials

@testset "CloudAtlas Tests" begin
    @time @testset "Polynomials" begin
        @test polyparity(Polynomial([1])) == 1
    end

    @time @testset "Fourier Modes" begin
        @test_nowarn FourierMode()
    end

    @time @testset "Basis Components" begin
        @test_nowarn BasisComponent()
        f = FourierMode(1, 1, 1)
        b = BasisComponent(1, f, f, Polynomial([1, 0, -1], :y), 0)
        @test Polynomial([1]) * b == b
    end
end
