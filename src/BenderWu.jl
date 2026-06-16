module BenderWu

# Bender-Wu method for perturbative energy levels of polynomial potentials
# Reference: arXiv:1608.08256

export Potential, clear_cache!
export max_k, A_kl, ε_l
export initialize_Akl_eps, fill_Akl!, eigenstate_coeffs
export find_epoly, epoly_taylor_derivatives, evaluate_epoly

"""
    Potential(vcoeffs)
    Potential(pairs)

Represents a polynomial potential with coefficients `vcoeffs`, where `vcoeffs[n]`
is the coefficient of x^(n+1). The index-to-power map is therefore offset by one:

| index     | `vcoeffs[1]` | `vcoeffs[2]` | `vcoeffs[3]` | … | `vcoeffs[n]` |
|-----------|--------------|--------------|--------------|---|--------------|
| power     | x²           | x³           | x⁴           | … | x^(n+1)      |

A second constructor accepts `power => coefficient` pairs for cases where the
positional convention is awkward — e.g. potentials with widely separated terms.
`Potential([2 => 0.5, 4 => 1.0])` is equivalent to `Potential([0.5, 0.0, 1.0])`.
Powers must be ≥ 2; duplicate powers are summed.

Carries its own memoization caches, which are GC-managed — create one instance
per potential and reuse it across calls.

`Rational{Int64}` (and any `Rational{<:Base.BitInteger}`) coefficients are
automatically promoted to `Rational{BigInt}` to prevent integer overflow at
higher perturbation orders.

The `vcoeffs` field must not be mutated after construction: caches are keyed
on these coefficients, and in-place changes would silently invalidate every
cached value. Construct a new `Potential` for a different polynomial.

Fields prefixed with an underscore (`_Akl_cache`, `_εl_cache`, `_cache_lock`)
are internal implementation details. They are not part of the public API and
may change without notice.

Cache access is guarded by a `ReentrantLock`, so a single `Potential` may be
shared safely across threads.

# Example
```julia
pot   = Potential([0.5, 0.0, 1.0])        # Float64, V(x) = 0.5x² + x⁴
pot   = Potential([2 => 0.5, 4 => 1.0])   # same potential via power => coeff pairs
pot_r = Potential([1//2, 0//1, 1//1])     # Rational — auto-promoted to Rational{BigInt}
epoly = find_epoly(2, pot)
```
"""
struct Potential{T}
    vcoeffs::Vector{T}
    ω::T
    _Akl_cache::Dict{Tuple{Int,Int,Int}, T}
    _εl_cache::Dict{Tuple{Int,Int}, T}
    _cache_lock::ReentrantLock
end

function Potential(vcoeffs::AbstractVector{T}) where T
    # The algorithm's domain is real-valued coefficients (ω = √(2·vcoeffs[1])
    # must be real). Catch non-real eltypes here so callers get a clear error
    # rather than a downstream MethodError from `>` on Complex/etc.
    T <: Real || throw(ArgumentError(
        "Potential supports real-valued coefficient types only " *
        "(Float64, BigFloat, Rational); got eltype = $T"))
    isempty(vcoeffs) && throw(ArgumentError("vcoeffs must be non-empty"))
    # Bender-Wu requires a positive harmonic baseline ω = √(2·vcoeffs[1]).
    # vcoeffs[1] ≤ 0 either makes ω complex (negative) or zero (which would
    # cause division-by-zero in the recursion below).
    vcoeffs[1] > zero(T) || throw(ArgumentError(
        "vcoeffs[1] (coefficient of x²) must be strictly positive; got $(vcoeffs[1])"))
    # Reuse the input directly when it is already a concrete Vector{T}; only
    # copy/convert for AbstractVector inputs (views, ranges, mismatched eltype).
    # The "must not be mutated after construction" contract in the docstring
    # covers cache integrity on the caller side.
    vc = vcoeffs isa Vector{T} ? vcoeffs : collect(T, vcoeffs)
    Potential(
        vc,
        _compute_ω(first(vcoeffs)),
        Dict{Tuple{Int,Int,Int}, T}(),
        Dict{Tuple{Int,Int}, T}(),
        ReentrantLock(),
    )
end

# Promote fixed-width rational coefficients to Rational{BigInt} to prevent
# integer overflow at higher perturbation orders.
Potential(vcoeffs::AbstractVector{Rational{T}}) where {T <: Base.BitInteger} =
    Potential(Rational{BigInt}.(vcoeffs))

function Potential(pairs::AbstractVector{<:Pair{<:Integer,T}}) where T
    isempty(pairs) && throw(ArgumentError("pairs must be non-empty"))
    any(first(p) < 2 for p in pairs) &&
        throw(ArgumentError("powers must be ≥ 2 (no x⁰ or x¹ terms allowed)"))
    maxp = maximum(first, pairs)
    vc = zeros(T, maxp - 1)
    for (p, c) in pairs
        vc[p - 1] += c
    end
    return Potential(vc)
end

"""
    clear_cache!(pot)

Empty the memoization caches inside `pot` and return `pot`. Useful when
sweeping over many `(ν, l)` values and you want to bound resident memory
between batches.
"""
function clear_cache!(pot::Potential)
    lock(pot._cache_lock) do
        empty!(pot._Akl_cache)
        empty!(pot._εl_cache)
    end
    return pot
end

"""
    max_k(pot, ν, l)

Return the maximum k-index K_l^(ν) at perturbation order `l` for quantum number `ν`.

K_l^(ν) bounds the support of the wave function coefficients A_{k,l}^(ν): all
coefficients with k > K_l^(ν) vanish. This bound depends on the degree of the
leading perturbation term in `pot`.
"""
function max_k(pot::Potential, ν::Int, l::Int)
    vcoeffs = pot.vcoeffs
    L = findfirst(!iszero, @view vcoeffs[2:end])
    # Pure harmonic oscillator: no perturbation terms, all higher-order
    # corrections vanish, so Kl = ν for l=0 and 0 otherwise.
    isnothing(L) && return iszero(l) ? ν : 0
    l < L ? (iszero(l) ? ν : 0) : ν + (L + 2) * (l ÷ L) + l % L
end

function _compute_ω(v::Rational{T}) where T
    two_v = 2 * v
    n, d = numerator(two_v), denominator(two_v)
    sn, sd = isqrt(n), isqrt(d)
    sn^2 == n && sd^2 == d ||
        throw(ArgumentError("2·vcoeffs[1] = $two_v is not a perfect rational square; ω is irrational"))
    return Rational{T}(sn, sd)
end
_compute_ω(v) = sqrt(2 * v)

# Shared update kernels used by both the recursive (A_kl/ε_l) and iterative
# (fill_Akl!) paths. They take accessor callables so the same formula can read
# values from a cache or from a pre-filled array. The compiler specialises on
# the callable types (FA, FE), so the closures inline without dynamic dispatch.

@inline function _akl_update(getA::FA, getε::FE, vcoeffs::Vector{T}, ω::T,
                              ν::Int, k::Int, l::Int) where {FA, FE, T}
    Akl = (k+2) * (k+1) * getA(k+2, l)
    # When k > ν and l > 0, the n == l term would call ε_l(ν, l), which itself
    # depends on A(ν+2, l) — still being computed in the recursive path. The
    # term vanishes anyway since A(k, 0) = 0 for k > ν, so omitting it is both
    # correct and circular-call-safe. (Same logic the iterative path expresses
    # as `for n = 1:l-1` in Step 1.)
    n_eps_max = (k > ν && l > 0) ? l - 1 : l
    for n = 1:n_eps_max
        Akl += 2 * getε(n) * getA(k, l-n)
    end
    for n = 1:l
        n+1 > length(vcoeffs) && continue
        iszero(vcoeffs[n+1]) && continue
        Akl += -2 * vcoeffs[n+1] * getA(k-n-2, l-n)
    end
    return Akl / (2 * ω * (k - ν))
end

@inline function _eps_update(getA::FA, vcoeffs::Vector{T},
                              ν::Int, l::Int) where {FA, T}
    ε = -(ν+2) * (ν+1) ÷ 2 * getA(ν+2, l)
    for n = 1:l
        n+1 > length(vcoeffs) && continue
        iszero(vcoeffs[n+1]) && continue
        ε += vcoeffs[n+1] * getA(ν-n-2, l-n)
    end
    return ε
end

"""
    A_kl(pot, ν, k, l)

Return the wave function expansion coefficient A_{k,l}^(ν).

These are the coefficients of the perturbative expansion of the ν-th eigenstate
in the harmonic oscillator basis, at perturbation order `l`. They satisfy a
recursive relation coupling different orders and indices; results are cached
inside `pot`.

Boundary conditions:
- Returns `zero` for k < 0, l < 0, or k > K_l^(ν)
- Returns `one` for k == ν, l == 0 (normalisation)
"""
function A_kl(pot::Potential{T}, ν::Int, k::Int, l::Int) where T
    # Cheap boundary cases — not worth caching
    if k < 0 || l < 0 return zero(T) end
    if k > ν && iszero(l) return zero(T) end
    if k == ν && iszero(l) return one(T) end
    if k == ν && l > 0 return zero(T) end
    if k > max_k(pot, ν, l) return zero(T) end

    lock(pot._cache_lock) do
        get!(pot._Akl_cache, (ν, k, l)) do
            _akl_update(
                (kk, ll) -> A_kl(pot, ν, kk, ll),
                n        -> ε_l(pot, ν, n),
                pot.vcoeffs, pot.ω, ν, k, l,
            )
        end
    end
end

"""
    ε_l(pot, ν, l)

Return the perturbative energy correction ε_l^(ν) for quantum number `ν` at
perturbation order `l`.

Odd orders vanish identically. Order l=0 returns the unperturbed harmonic
energy ω·(ν + 1/2), where ω = √(2·vcoeffs[1]). Results for even orders are
cached inside `pot`.

# `l` is the BW recursion index, not the physical perturbation order

Let `L` be the index of the first non-zero entry in `vcoeffs[2:end]` (so `L = 2`
for a quartic potential, `L = 4` for a sextic, `L = 6` for an octic). Bender–Wu
indexes corrections by the recursion order `l`, related to the physical
perturbation order `k` by `l = k·L`. Consequently `ε_l(pot, ν, l) = 0` whenever
`l mod L ≠ 0` — for a sextic, `ε_l(pot, ν, 2) = 0` and the first non-zero
correction sits at `l = 4`. See the "Energy corrections" section of the README
for a worked sextic example.
"""
function ε_l(pot::Potential{T}, ν::Int, l::Int) where T
    # Cheap boundary cases — not worth caching
    if isodd(l) return zero(T) end
    ω = pot.ω
    if iszero(l) return ω * (ν + one(T)/2) end

    lock(pot._cache_lock) do
        get!(pot._εl_cache, (ν, l)) do
            _eps_update((kk, ll) -> A_kl(pot, ν, kk, ll), pot.vcoeffs, ν, l)
        end
    end
end

"""
    initialize_Akl_eps(pot, ν, l)

Allocate and return zero-initialised arrays `(Akl, ε)` sized for the iterative
computation up to perturbation order `l` for quantum number `ν`.

`Akl` has dimensions `(K_l^(ν) + 3) × (l + 1)` and `ε` has length `l + 1`.
Element type matches `eltype(pot.vcoeffs)`. Pass these arrays to `fill_Akl!`.
"""
function initialize_Akl_eps(pot::Potential, ν::Int, l::Int)
    # max(ν, …): pure harmonic potentials have max_k = 0 for l > 0, but the
    # boundary condition still writes Akl[ν+1, 1].
    kmax = max(ν, max_k(pot, ν, l))
    T = eltype(pot.vcoeffs)
    return zeros(T, kmax+3, l+1), zeros(T, l+1)
end

"""
    fill_Akl!(Akl, ε, pot, ν, maxorder)

Fill pre-allocated arrays `Akl` and `ε` in-place with wave function coefficients
and energy corrections up to perturbation order `maxorder` for quantum number `ν`.

This is a non-recursive, type-stable alternative to the cached `A_kl`/`ε_l`
functions. Each order `l` is computed in three steps:

1. Compute `Akl[k, l]` for k > ν (descending from K_l^(ν))
2. Compute `ε[l]` from the boundary condition at k = ν
3. Compute `Akl[k, l]` for k < ν (descending from ν−1)

Use `initialize_Akl_eps` to allocate arrays of the correct size.

# Note
Array indexing is 1-based: `Akl[k+1, l+1]` holds the coefficient for index k
at order l, and `ε[l+1]` holds the energy correction at order l.
"""
function fill_Akl!(Akl, ε, pot::Potential, ν::Int, maxorder::Int)
    vcoeffs = pot.vcoeffs
    ω = pot.ω
    T = eltype(vcoeffs)
    Akl[ν+1, 1] = one(ω)
    ε[1] = ω * (ν + one(T)/2)
    # Accessors give _akl_update / _eps_update read-only views of Akl/ε that
    # match the boundary semantics of the recursive path (k < 0 ⇒ zero).
    # Array indexing is 1-based: A(k, l) lives at Akl[k+1, l+1], ε(n) at ε[n+1].
    getA = (k, l) -> k < 0 ? zero(T) : Akl[k+1, l+1]
    getε = n -> ε[n+1]
    for l = 0:maxorder
        kmax = max_k(pot, ν, l)
        if l > 0
            for k = kmax:-1:ν+1
                Akl[k+1, l+1] = _akl_update(getA, getε, vcoeffs, ω, ν, k, l)
            end
            ε[l+1] = _eps_update(getA, vcoeffs, ν, l)
        end
        for k = ν-1:-1:0
            Akl[k+1, l+1] = _akl_update(getA, getε, vcoeffs, ω, ν, k, l)
        end
    end
    nothing
end

"""
    eigenstate_coeffs(pot, ν, l)

Return the 2D array of wave function expansion coefficients A_{k,l}^(ν) up to
BW recursion order `l` for quantum number `ν`.

The ν-th perturbed eigenstate is expanded as ψ_ν(x) = e^{-ω x² / 2}·∑_l λ^l ∑_k
A_{k,l}^(ν) · x^k (λ is the coupling). This function returns the full coefficient
table A_{k,l}^(ν), packed into a single 2D array indexed as `Akl[k+1, l+1]`
(1-based offset matches [`fill_Akl!`](@ref)). The k=ν entry of the zeroth-order
column is fixed by normalisation (`Akl[ν+1, 1] == one(T)`); other entries in
that column encode the Hermite-like polynomial of the unperturbed eigenstate
and are generally non-zero for ν ≥ 2.

This is a convenience wrapper around [`initialize_Akl_eps`](@ref) followed by
[`fill_Akl!`](@ref) for callers who want only the eigenstate coefficients and
not the energy-correction array. Element type matches `eltype(pot.vcoeffs)`.

# `l` is the BW recursion index
See the [`ε_l`](@ref) docstring for the relation `l = k·L` between the BW
recursion index and the physical perturbation order.

# Example
```julia
pot = Potential([0.5, 0.0, 1.0])      # V = x²/2 + x⁴
Akl = eigenstate_coeffs(pot, 0, 4)    # ground state, up to BW order 4
Akl[1, 1] == 1.0                       # unperturbed ground state coefficient
```
"""
function eigenstate_coeffs(pot::Potential, ν::Int, l::Int)
    Akl, ε = initialize_Akl_eps(pot, ν, l)
    fill_Akl!(Akl, ε, pot, ν, l)
    return Akl
end

"""
    find_epoly(order, pot)

Return the coefficients of the energy polynomial ε^(order)(ν) at perturbation
order `order`.

The energy eigenvalue at perturbation order `order` is a polynomial in the
quantum number ν. This function evaluates ε_l^(ν) at `order/2 + 2` integer
values of ν, computes divided differences (the Newton form), and converts to
the monomial basis. Compared to a Vandermonde solve this is O(n²) instead of
O(n³), exact in rational arithmetic, and numerically stable for `Float64` at
high orders.

Returns a zero vector for odd `order` (all odd-order corrections vanish).
The element type matches `eltype(pot.vcoeffs)`, so pass a `BigFloat`-based
`Potential` for arbitrary-precision results.

# Example
```julia
pot   = Potential([0.5, 0.0, 1.0])
epoly = find_epoly(2, pot)          # first correction for quartic oscillator
evaluate_epoly(3, epoly)            # energy correction at ν = 3
```
"""
function find_epoly(order::Int, pot::Potential)
    T = eltype(pot.vcoeffs)
    if isodd(order)
        # Size matches the even-order layout (order ÷ 2 + 2) so callers can
        # index uniformly across orders; entries are zero since odd-order
        # corrections vanish identically.
        return zeros(T, order ÷ 2 + 2)
    end
    n = order ÷ 2 + 2
    dd = [ε_l(pot, k, order) for k = 0:n-1]
    # In-place divided differences for unit-spaced nodes 0..n-1: divide by k.
    for k = 1:n-1, i = n:-1:k+1
        dd[i] = (dd[i] - dd[i-1]) / k
    end
    # Convert Newton form ∑ dd[k+1] · x(x-1)…(x-k+1) to monomial coefficients.
    # In-place Horner-like update: walking j from len downward writes
    # coeffs[j+1] += coeffs[j] using still-fresh coeffs[j], then overwrites
    # coeffs[j] = -(k-1)*coeffs[j]. The next step (smaller j) supplies the
    # additive term for the now-overwritten coeffs[j].
    coeffs = zeros(T, n)
    coeffs[1] = dd[n]
    len = 1
    for k = n-1:-1:1
        for j = len:-1:1
            coeffs[j+1] += coeffs[j]
            coeffs[j]    = -(k-1) * coeffs[j]
        end
        coeffs[1] += dd[k]
        len += 1
    end
    return coeffs
end

"""
    epoly_taylor_derivatives(epoly)

Return the derivatives of the energy polynomial `epoly` evaluated at ν = 0.

Given `epoly = [c_0, c_1, …, c_n]` representing the polynomial ∑ c_k · ν^k,
the k-th output element (1-indexed) is k! · c_k, i.e. the k-th derivative of
the polynomial at ν = 0. Output length is `length(epoly) - 1` (the constant
term c_0 is not a derivative).
"""
function epoly_taylor_derivatives(epoly)
    isempty(epoly) && throw(ArgumentError("epoly must be non-empty"))
    T = typeof(epoly[1])
    ds = Array{T}(undef, length(epoly)-1)
    for k in 1:length(epoly)-1
        # factorial(20) fits in Int64 but factorial(21) overflows — switch to
        # BigInt at the threshold to avoid OverflowError while keeping the
        # cheap path for typical (small-order) inputs.
        fact_k = k < 20 ? factorial(k) : factorial(big(k))
        ds[k] = fact_k * epoly[k+1]
    end
    return ds
end

"""
    evaluate_epoly(n, epoly)

Evaluate the energy polynomial `epoly` at quantum number `n`.

`epoly` is a coefficient vector [c_0, c_1, ..., c_m] representing
∑_k c_k · n^(k-1).
"""
function evaluate_epoly(n::Int, epoly)
    isempty(epoly) && throw(ArgumentError("epoly must be non-empty"))
    T = typeof(epoly[1])
    res = zero(T)
    for k in length(epoly):-1:1
        res = res * n + epoly[k]
    end
    return res
end

end # module BenderWu
