using Test
using BenderWu

@testset "BenderWu" begin

    pot = Potential([0.5, 0.0, 1.0])  # quartic oscillator: V = x²/2 + x⁴, ω = 1

    @testset "Zeroth-order energy (harmonic)" begin
        # ε_l(pot, ν, 0) = ω*(ν + 1/2) = ν + 0.5
        @test ε_l(pot, 0, 0) ≈ 0.5
        @test ε_l(pot, 1, 0) ≈ 1.5
        @test ε_l(pot, 3, 0) ≈ 3.5
    end

    @testset "Odd perturbation orders vanish" begin
        @test iszero(ε_l(pot, 0, 1))
        @test iszero(ε_l(pot, 2, 3))
    end

    @testset "First perturbative correction (order l=2)" begin
        # Analytical result: ⟨ν|x⁴|ν⟩ = (6ν² + 6ν + 3) / 4
        @test ε_l(pot, 0, 2) ≈ 3/4
        @test ε_l(pot, 1, 2) ≈ 15/4
        @test ε_l(pot, 2, 2) ≈ 39/4
    end

    @testset "evaluate_epoly round-trips find_epoly" begin
        epoly = find_epoly(2, pot)
        @test evaluate_epoly(0, epoly) ≈ ε_l(pot, 0, 2)
        @test evaluate_epoly(1, epoly) ≈ ε_l(pot, 1, 2)
        @test evaluate_epoly(2, epoly) ≈ ε_l(pot, 2, 2)
    end

    @testset "Float64 and BigFloat caches are independent" begin
        pot_bf = Potential(BigFloat.([0.5, 0.0, 1.0]))
        @test ε_l(pot_bf, 0, 2) ≈ big(3)/4
        @test ε_l(pot_bf, 1, 2) ≈ big(15)/4
        # Ensure the Float64 cache is unaffected
        @test ε_l(pot, 0, 2) ≈ 3/4
    end

    @testset "Rational arithmetic (exact results)" begin
        pot_r = Potential([1//2, 0//1, 1//1])
        @test ε_l(pot_r, 0, 0) == 1//2
        @test ε_l(pot_r, 1, 0) == 3//2
        @test ε_l(pot_r, 0, 2) == 3//4
        @test ε_l(pot_r, 1, 2) == 15//4
        @test ε_l(pot_r, 2, 2) == 39//4
    end

end
