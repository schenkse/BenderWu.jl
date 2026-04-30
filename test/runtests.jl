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
        # Output element type must match the polynomial's element type
        @test evaluate_epoly(2, epoly) isa Float64
        epoly_bf = find_epoly(2, Potential(BigFloat.([0.5, 0.0, 1.0])))
        @test evaluate_epoly(2, epoly_bf) isa BigFloat
        epoly_r = find_epoly(2, Potential([1//2, 0//1, 1//1]))
        @test evaluate_epoly(2, epoly_r) isa Rational{BigInt}
    end

    @testset "Float64 and BigFloat caches are independent" begin
        pot_bf = Potential(BigFloat.([0.5, 0.0, 1.0]))
        @test ε_l(pot_bf, 0, 2) ≈ big(3)/4
        @test ε_l(pot_bf, 1, 2) ≈ big(15)/4
        # Ensure the Float64 cache is unaffected
        @test ε_l(pot, 0, 2) ≈ 3/4
    end

    @testset "fill_Akl! matches recursive ε_l" begin
        # Iterative and recursive implementations must agree on energy corrections
        for ν in 0:3
            maxorder = 4
            Akl, ε_arr = initialize_Akl_eps(pot, ν, maxorder)
            fill_Akl!(Akl, ε_arr, pot, ν, maxorder)
            for l in 0:maxorder
                @test ε_arr[l+1] ≈ ε_l(pot, ν, l)
            end
        end
    end

    @testset "epoly_taylor_derivatives" begin
        # For quartic V = x²/2 + x⁴, the order-2 energy polynomial is
        # ε^(2)(ν) = (6ν² + 6ν + 3)/4 = (3/2)ν² + (3/2)ν + 3/4
        # so epoly = [3/4, 3/2, 3/2] and the derivatives at ν=0 are:
        # ds[1] = 1! * (3/2) = 1.5  (first derivative at ν=0)
        # ds[2] = 2! * (3/2) = 3.0  (second derivative at ν=0)
        epoly = find_epoly(2, pot)
        ds = epoly_taylor_derivatives(epoly)
        @test length(ds) == length(epoly) - 1
        @test ds[1] ≈ 1.5
        @test ds[2] ≈ 3.0
    end

    @testset "Pure harmonic potential" begin
        # No perturbation terms: all higher-order corrections vanish, only
        # ε_l(pot, ν, 0) = ω*(ν + 1/2) is non-zero.
        pot_h = Potential([0.5])  # V = x²/2, ω = 1
        for ν in 0:5
            @test ε_l(pot_h, ν, 0) ≈ ν + 0.5
            @test iszero(ε_l(pot_h, ν, 2))
            @test iszero(ε_l(pot_h, ν, 4))
        end
        # fill_Akl! must not go out of bounds for ν ≥ 3 where max_k(·, l>0) = 0
        for ν in 0:5
            maxorder = 4
            Akl, ε_arr = initialize_Akl_eps(pot_h, ν, maxorder)
            fill_Akl!(Akl, ε_arr, pot_h, ν, maxorder)
            @test ε_arr[1] ≈ ν + 0.5
            for l in 1:maxorder
                @test iszero(ε_arr[l+1])
            end
        end
    end

    @testset "Rational arithmetic (exact results)" begin
        pot_r = Potential([1//2, 0//1, 1//1])
        # Integer-typed rationals are auto-promoted to Rational{BigInt}
        @test eltype(pot_r.vcoeffs) == Rational{BigInt}
        @test ε_l(pot_r, 0, 0) == 1//2
        @test ε_l(pot_r, 1, 0) == 3//2
        @test ε_l(pot_r, 0, 2) == 3//4
        @test ε_l(pot_r, 1, 2) == 15//4
        @test ε_l(pot_r, 2, 2) == 39//4
        # Order 14 must complete without OverflowError
        @test length(find_epoly(14, pot_r)) == 9
    end

end
