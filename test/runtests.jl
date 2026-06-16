using Test
using BenderWu
using Base.Threads

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

    @testset "Empty epoly throws ArgumentError" begin
        # Without the guard, both fall through to `typeof(epoly[1])` and emit
        # a confusing BoundsError. Empty input is invalid; surface it clearly.
        @test_throws ArgumentError evaluate_epoly(0, Float64[])
        @test_throws ArgumentError epoly_taylor_derivatives(Float64[])
    end

    @testset "Constructor input validation" begin
        @test_throws ArgumentError Potential(Float64[])
        @test_throws ArgumentError Potential([0.0, 0.0, 1.0])   # no x² term
        @test_throws ArgumentError Potential([-0.5, 0.0, 1.0])  # negative ω²
        @test_throws ArgumentError Potential([0//1, 0//1, 1//1])
        # Non-real eltypes are rejected upfront with a clear ArgumentError
        # rather than a downstream MethodError on `>`.
        @test_throws ArgumentError Potential(ComplexF64[1.0 + 0im, 0.0, 1.0])
    end

    @testset "Vector{T} input is reused, AbstractVector input is copied" begin
        # Common case: passing a Vector{T} of the right eltype reuses the
        # buffer (no defensive copy).
        src = [0.5, 0.0, 1.0]
        @test Potential(src).vcoeffs === src
        # AbstractVector inputs (views, ranges) still go through collect.
        v = @view src[1:end]
        pot_v = Potential(v)
        @test pot_v.vcoeffs !== src
        @test pot_v.vcoeffs == src
        @test pot_v.vcoeffs isa Vector{Float64}
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
        epoly14 = find_epoly(14, pot_r)
        @test length(epoly14) == 9
        # Ground-state (ν=0) corrections of V = x²/2 + x⁴ are tabulated in
        # Bender & Wu, Phys. Rev. D 7, 1620 (1973). With BW recursion index
        # l = k·L (L=2 for quartic), l=14 is physical order 7.
        # Pin several literature coefficients via the constant term of
        # find_epoly(l, pot_r), which equals ε_l(pot_r, 0, l).
        @test find_epoly(8,  pot_r)[1] == -30885 // 128
        @test find_epoly(10, pot_r)[1] == 916731 // 256
        @test find_epoly(12, pot_r)[1] == -65518401 // 1024
        @test epoly14[1]               == 2723294673 // 2048
    end

    @testset "A_kl boundary values" begin
        # Documented boundary conditions on A_{k,l}^(ν).
        for ν in 0:3
            # Normalisation at order 0.
            @test A_kl(pot, ν, ν, 0) == 1.0
            # k > ν vanishes at order 0.
            for k in ν+1:ν+3
                @test A_kl(pot, ν, k, 0) == 0.0
            end
            # k = ν vanishes at every higher order.
            for l in 1:4
                @test A_kl(pot, ν, ν, l) == 0.0
            end
            # k < 0 and l < 0 vanish.
            @test A_kl(pot, ν, -1, 0) == 0.0
            @test A_kl(pot, ν, 0, -1) == 0.0
            # k > K_l^(ν) vanishes.
            @test A_kl(pot, ν, max_k(pot, ν, 2) + 1, 2) == 0.0
        end
    end

    @testset "BW recursion index l=2 vanishes for sextic/octic" begin
        # l is the BW recursion index, not the physical perturbation order.
        # For leading perturbation degree L, ε_l = 0 whenever l mod L ≠ 0.
        # Sextic: L=4 — l=2 vanishes, l=4 is the first physical correction.
        # Octic:  L=6 — l=2 and l=4 vanish, l=6 is the first correction.
        pot6 = Potential([0.5, 0.0, 0.0, 0.0, 1.0])
        pot8 = Potential([0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
        for ν in 0:3
            @test iszero(ε_l(pot6, ν, 2))
            @test iszero(ε_l(pot8, ν, 2))
            @test iszero(ε_l(pot8, ν, 4))
        end
        @test !iszero(ε_l(pot6, 0, 4))
        @test !iszero(ε_l(pot8, 0, 6))
    end

    @testset "find_epoly agrees across BigFloat and Rational precision" begin
        # Order-10 quartic energy polynomial must match to ~30 decimal digits
        # whether computed in BigFloat or exact rational arithmetic.
        pot_bf = Potential(BigFloat.([0.5, 0.0, 1.0]))
        pot_r  = Potential([1//2, 0//1, 1//1])
        e_bf   = find_epoly(10, pot_bf)
        e_r    = find_epoly(10, pot_r)
        @test length(e_bf) == length(e_r)
        for (a, b) in zip(e_bf, e_r)
            @test isapprox(a, BigFloat(b); rtol = BigFloat(10)^-30)
        end
    end

    @testset "clear_cache! empties the caches" begin
        pot_x = Potential([0.5, 0.0, 1.0])
        ε_l(pot_x, 0, 4)
        A_kl(pot_x, 1, 3, 2)
        @test !isempty(pot_x._εl_cache)
        @test !isempty(pot_x._Akl_cache)
        @test clear_cache!(pot_x) === pot_x
        @test isempty(pot_x._εl_cache)
        @test isempty(pot_x._Akl_cache)
        # Results must still be reproducible after clearing.
        @test ε_l(pot_x, 0, 4) ≈ ε_l(Potential([0.5, 0.0, 1.0]), 0, 4)
    end

    @testset "ε_l memoization: repeated calls hit the cache" begin
        # Reaches into the private _εl_cache to assert the memoization
        # contract: calling ε_l with identical arguments must not grow the
        # cache after the first call. Tied to A2's lock-guarded get!.
        pot_c = Potential([0.5, 0.0, 1.0])
        ε_l(pot_c, 0, 4)
        n_before = length(pot_c._εl_cache)
        for _ in 1:5
            ε_l(pot_c, 0, 4)
        end
        @test length(pot_c._εl_cache) == n_before
    end

    @testset "Sextic potential (V = x²/2 + x⁶)" begin
        pot6 = Potential([0.5, 0.0, 0.0, 0.0, 1.0])
        # find_epoly + evaluate_epoly must round-trip ε_l for several (ν, l).
        for l in (0, 2, 4)
            epoly = find_epoly(l, pot6)
            for ν in 0:3
                @test evaluate_epoly(ν, epoly) ≈ ε_l(pot6, ν, l)
            end
        end
        # Iterative path must match recursive path.
        maxorder = 4
        for ν in 0:3
            Akl, ε_arr = initialize_Akl_eps(pot6, ν, maxorder)
            fill_Akl!(Akl, ε_arr, pot6, ν, maxorder)
            for l in 0:maxorder
                @test ε_arr[l+1] ≈ ε_l(pot6, ν, l)
            end
        end
    end

    @testset "Octic potential (V = x²/2 + x⁸)" begin
        pot8 = Potential([0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
        for l in (0, 2, 4)
            epoly = find_epoly(l, pot8)
            for ν in 0:2
                @test evaluate_epoly(ν, epoly) ≈ ε_l(pot8, ν, l)
            end
        end
    end

    @testset "Mixed-parity potential (V = x²/2 + x³ + x⁴)" begin
        # Highest power must be even (physics); odd-power *terms* are allowed.
        # Leading non-zero perturbation index is L=1, exercising a different
        # branch of max_k than the quartic (where L=2).
        pot_m = Potential([0.5, 1.0, 1.0])
        maxorder = 4
        for ν in 0:3
            Akl, ε_arr = initialize_Akl_eps(pot_m, ν, maxorder)
            fill_Akl!(Akl, ε_arr, pot_m, ν, maxorder)
            for l in 0:maxorder
                @test ε_arr[l+1] ≈ ε_l(pot_m, ν, l)
            end
        end
        # Odd orders still vanish identically.
        @test iszero(ε_l(pot_m, 0, 1))
        @test iszero(ε_l(pot_m, 2, 3))
    end

    @testset "Thread-safe cache (shared Potential, many threads)" begin
        # Race-detection smoke test: hammer ε_l from many threads on a shared
        # Potential. Without the ReentrantLock this corrupts the get! Dict and
        # surfaces as wrong values or rare exceptions.
        pot_t = Potential([0.5, 0.0, 1.0])
        expected = [ε_l(Potential([0.5, 0.0, 1.0]), ν, 2*l) for ν in 0:3, l in 0:5]
        results = Array{Float64}(undef, 4, 6, 32)
        @threads for trial in 1:32
            for ν in 0:3, l in 0:5
                results[ν+1, l+1, trial] = ε_l(pot_t, ν, 2*l)
            end
        end
        for trial in 1:32, ν in 0:3, l in 0:5
            @test results[ν+1, l+1, trial] == expected[ν+1, l+1]
        end
    end

    @testset "find_epoly at higher order (BigFloat)" begin
        pot_bf = Potential(BigFloat.([0.5, 0.0, 1.0]))
        epoly = find_epoly(20, pot_bf)
        @test eltype(epoly) == BigFloat
        for ν in 0:3
            @test evaluate_epoly(ν, epoly) ≈ ε_l(pot_bf, ν, 20)
        end
    end

    @testset "Potential(pairs) constructor" begin
        # Equivalent to the quartic Potential([0.5, 0.0, 1.0]).
        pot_p = Potential([2 => 0.5, 4 => 1.0])
        @test pot_p.vcoeffs == [0.5, 0.0, 1.0]
        @test ε_l(pot_p, 0, 2) ≈ 3/4
        @test ε_l(pot_p, 1, 2) ≈ 15/4
        # Order of pairs and gaps in powers should not matter.
        @test Potential([4 => 1.0, 2 => 0.5]).vcoeffs == [0.5, 0.0, 1.0]
        # Duplicate powers accumulate.
        @test Potential([2 => 0.25, 2 => 0.25, 4 => 1.0]).vcoeffs == [0.5, 0.0, 1.0]
        # Validation errors.
        @test_throws ArgumentError Potential(Pair{Int,Float64}[])
        @test_throws ArgumentError Potential([1 => 1.0, 2 => 0.5])  # x¹ term
        @test_throws ArgumentError Potential([2 => -0.5, 4 => 1.0]) # bad ω
    end

    @testset "Iterative path with Potential(pairs) constructor" begin
        # fill_Akl! must agree with recursive ε_l when the Potential was built
        # from the pairs constructor (covers the iterative path on a non-Vector
        # construction route).
        pot_p = Potential([2 => 0.5, 4 => 1.0])
        maxorder = 6
        for ν in 0:3
            Akl, ε_arr = initialize_Akl_eps(pot_p, ν, maxorder)
            fill_Akl!(Akl, ε_arr, pot_p, ν, maxorder)
            for l in 0:maxorder
                @test ε_arr[l+1] ≈ ε_l(pot_p, ν, l)
            end
        end
    end

    @testset "eigenstate_coeffs returns the A_kl matrix" begin
        # Thin wrapper over initialize_Akl_eps + fill_Akl!. Must match the
        # iterative-path layout exactly: Akl[k+1, l+1] = A_kl(pot, ν, k, l).
        maxorder = 4
        for ν in 0:3
            Akl = eigenstate_coeffs(pot, ν, maxorder)
            # Size matches initialize_Akl_eps.
            kmax = max(ν, max_k(pot, ν, maxorder))
            @test size(Akl) == (kmax + 3, maxorder + 1)
            # Normalisation: the k = ν entry of the zeroth-order column is 1.
            @test Akl[ν+1, 1] == 1.0
            # Full agreement with the recursive path on every (k, l) cell —
            # this also implicitly covers the Hermite-like structure of the
            # zeroth-order column for ν ≥ 2.
            for l in 0:maxorder, k in 0:kmax
                @test Akl[k+1, l+1] ≈ A_kl(pot, ν, k, l)
            end
        end

        # Element type follows eltype(pot.vcoeffs).
        @test eltype(eigenstate_coeffs(pot, 0, 2)) == Float64
        pot_bf = Potential(BigFloat.([0.5, 0.0, 1.0]))
        @test eltype(eigenstate_coeffs(pot_bf, 0, 2)) == BigFloat
        pot_r = Potential([1//2, 0//1, 1//1])
        Akl_r = eigenstate_coeffs(pot_r, 1, 4)
        @test eltype(Akl_r) == Rational{BigInt}
        # Exact agreement with the recursive path under rational arithmetic.
        for l in 0:4, k in 0:size(Akl_r, 1) - 1
            @test Akl_r[k+1, l+1] == A_kl(pot_r, 1, k, l)
        end

        # Sextic: first physical correction sits at BW recursion order l = 4.
        pot6 = Potential([0.5, 0.0, 0.0, 0.0, 1.0])
        Akl6 = eigenstate_coeffs(pot6, 0, 4)
        for l in 0:4, k in 0:size(Akl6, 1) - 1
            @test Akl6[k+1, l+1] ≈ A_kl(pot6, 0, k, l)
        end
    end

    @testset "Iterative path with Rational coefficients (exact)" begin
        # All prior iterative-path tests use Float64. Cover the iterative path
        # under exact Rational{BigInt} arithmetic: results must be type-stable
        # (no Float64 contamination) and exactly equal — not just ≈ — to the
        # recursive path.
        pot_r = Potential([1//2, 0//1, 1//1])
        maxorder = 6
        for ν in 0:3
            Akl, ε_arr = initialize_Akl_eps(pot_r, ν, maxorder)
            @test eltype(Akl) == Rational{BigInt}
            @test eltype(ε_arr) == Rational{BigInt}
            fill_Akl!(Akl, ε_arr, pot_r, ν, maxorder)
            for l in 0:maxorder
                @test ε_arr[l+1] == ε_l(pot_r, ν, l)
            end
        end
    end

    @testset "Irrational rational ω is rejected with ArgumentError" begin
        # 2·(1//3) = 2//3 is not a perfect rational square, so ω would be
        # irrational. Must surface as ArgumentError, consistent with every
        # other construction-time validation failure.
        @test_throws ArgumentError Potential([1//3, 0//1, 1//1])
    end

    @testset "find_epoly returns an all-zero vector at odd order" begin
        # Odd-order corrections vanish identically; find_epoly short-circuits
        # to zeros sized to match the next-lower even order (order ÷ 2 + 2).
        e3 = find_epoly(3, pot)
        @test length(e3) == 3            # == length(find_epoly(2, pot))
        @test all(iszero, e3)
        @test eltype(e3) == Float64
        # Exact arithmetic: odd order returns exact rational zeros.
        pot_r = Potential([1//2, 0//1, 1//1])
        e5 = find_epoly(5, pot_r)
        @test length(e5) == 4            # 5 ÷ 2 + 2
        @test all(iszero, e5)
        @test eltype(e5) == Rational{BigInt}
    end

    @testset "epoly_taylor_derivatives BigInt-factorial branch (k ≥ 20)" begin
        # Exercise the k ≥ 20 path (factorial(big(k))) without an expensive
        # high-order find_epoly: a synthetic length-22 epoly makes k run 1:21.
        # All-ones rational coefficients let us pin the exact factorial values.
        epoly_big = [big(1)//1 for _ in 1:22]
        ds = epoly_taylor_derivatives(epoly_big)
        @test length(ds) == 21
        @test ds[19] == factorial(19)        # last Int-factorial step  (k < 20)
        @test ds[20] == factorial(big(20))   # first BigInt step        (k ≥ 20)
        @test ds[21] == factorial(big(21))   # overflows Int64 — exact via BigInt
    end

end
