using Test
using BenderWu

@testset "BenderWu" begin

    vcoeffs = [0.5, 0.0, 1.0]  # quartic oscillator: V = x²/2 + x⁴, ω = 1

    @testset "Zeroth-order energy (harmonic)" begin
        # ε_l(ν, 0) = ω*(ν + 1/2) = ν + 0.5
        @test ε_l(0, 0, vcoeffs) ≈ 0.5
        @test ε_l(1, 0, vcoeffs) ≈ 1.5
        @test ε_l(3, 0, vcoeffs) ≈ 3.5
    end

    @testset "Odd perturbation orders vanish" begin
        @test iszero(ε_l(0, 1, vcoeffs))
        @test iszero(ε_l(2, 3, vcoeffs))
    end

    @testset "First perturbative correction (order l=2)" begin
        # Analytical result: ⟨ν|x⁴|ν⟩ = (6ν² + 6ν + 3) / 4
        @test ε_l(0, 2, vcoeffs) ≈ 3/4
        @test ε_l(1, 2, vcoeffs) ≈ 15/4
        @test ε_l(2, 2, vcoeffs) ≈ 39/4
    end

    @testset "evaluate_epoly round-trips find_epoly" begin
        epoly = find_epoly(2, vcoeffs)
        @test evaluate_epoly(0, epoly) ≈ ε_l(0, 2, vcoeffs)
        @test evaluate_epoly(1, epoly) ≈ ε_l(1, 2, vcoeffs)
        @test evaluate_epoly(2, epoly) ≈ ε_l(2, 2, vcoeffs)
    end

end
