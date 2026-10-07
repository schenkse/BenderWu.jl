# BenderWu.jl

Julia implementation of the **Bender-Wu method** for computing perturbative energy corrections to eigenvalues of 1D quantum systems with polynomial potentials.

The Hamiltonian is H = -½∂ₓ² + V(x), where V(x) = Σ vcoeffs[n] · xⁿ⁺¹. The package expands the energy eigenvalue E_ν in a coupling `g` using the Bender-Wu recursive relations. See [Energy corrections](#energy-corrections-ε_lpot-ν-l) for how `g` enters H.

[![Julia ≥ 1.10](https://img.shields.io/badge/Julia-≥1.10-9558B2?logo=julia)](https://julialang.org)
[![No dependencies](https://img.shields.io/badge/dependencies-none-brightgreen)](Project.toml)

> **Built with LLMs:** This project was developed with the help of AI coding tools, primarily [Claude Code](https://claude.ai/code) by Anthropic. All code has been reviewed and is maintained by the author.

## References

- C. M. Bender & T. T. Wu, *Phys. Rev.* **184**, 1231 (1969) — original recursion relations for the anharmonic oscillator  
  <https://doi.org/10.1103/PhysRev.184.1231>
- C. M. Bender & T. T. Wu, *Phys. Rev. D* **7**, 1620 (1973) — extension to higher-order perturbation theory  
  <https://doi.org/10.1103/PhysRevD.7.1620>
- T. Sulejmanpasic & M. Ünsal, *arXiv:1608.08256* (2016) — algorithmic presentation with a Wolfram/Mathematica implementation; this Julia package follows those recursive relations  
  <https://arxiv.org/abs/1608.08256>

## Setup

### Install from GitHub

```julia
using Pkg
Pkg.add(url="https://github.com/schenkse/BenderWu.jl")
using BenderWu
```

### Clone for development

```bash
git clone https://github.com/schenkse/BenderWu.jl
cd BenderWu.jl
```

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using BenderWu
```

## Usage

### Creating a potential

A `Potential` is constructed from a coefficient vector where `vcoeffs[n]` is the coefficient of xⁿ⁺¹:

```julia
pot = Potential([0.5, 0.0, 1.0])   # V(x) = (1/2)x² + x⁴  (quartic oscillator, ω = 1)
```

The index-to-power offset is by one: `vcoeffs[1]` is the x² coefficient, `vcoeffs[2]` is the x³ coefficient, and so on. For potentials where this convention is awkward — e.g. with widely separated terms — pass `power => coefficient` pairs instead:

```julia
pot = Potential([2 => 0.5, 4 => 1.0])    # same V(x) as above
```

The frequency ω is derived automatically from the leading term: ω = √(2·vcoeffs[1]).

Create one `Potential` per potential and reuse it — results are memoized inside the struct and freed automatically when it goes out of scope. Cache access is guarded by a `ReentrantLock`, so a single `Potential` can be shared safely across threads. If you sweep many `(ν, l)` and want to bound resident memory between batches, call `clear_cache!(pot)`.

### Energy corrections ε_l(pot, ν, l)

The package expands in the coupling `g` of the scaled Hamiltonian

H(g) = -½∂ₓ² + V(gx)/g²,

where V is the polynomial from `vcoeffs`. A monomial `c·x^p` in V enters H(g) as `c·g^(p-2)·x^p`, so the harmonic term carries no `g`. `ε_l(pot, ν, l)` returns the coefficient of `g^l` in the energy E_ν(g) for quantum number `ν`.

```julia
ε_l(pot, 0, 0)   # → 0.5    (ground state, unperturbed: ω·(0 + 1/2))
ε_l(pot, 1, 0)   # → 1.5    (ν=1, unperturbed: ω·(1 + 1/2))
ε_l(pot, 0, 2)   # → 0.75   (first-order quartic correction ⟨0|x⁴|0⟩ = 3/4)
ε_l(pot, 1, 2)   # → 3.75
ε_l(pot, 2, 2)   # → 9.75
```

Odd orders vanish for every potential and are returned as exact zero. Flipping g → −g together with x → −x leaves H(g) unchanged, so E_ν(g) is even in `g`.

#### `l` is the BW recursion index, not the physical perturbation order

Which orders survive depends on the monomials in V. Call `p − 2` the index of an anharmonic monomial `x^p`.

**Single anharmonic monomial.** If V has one term `c·x^p` besides x², set `L = p − 2`. The physical coupling is `λ = g^L`, and the k-th order in λ sits at `l = k·L`. That gives `l = 2k` for a quartic, `l = 4k` for a sextic and `l = 6k` for an octic. `ε_l(pot, ν, l) = 0` unless `L` divides `l`.

**Worked sextic example** (V = x²/2 + x⁶, so `L = 4`):

```julia
pot6 = Potential([0.5, 0.0, 0.0, 0.0, 1.0])
ε_l(pot6, 0, 2)   # → 0.0     (no physical correction here — l = 2 is below L)
ε_l(pot6, 0, 4)   # ≠ 0       (first physical correction)
```

**Mixed potentials.** With several anharmonic monomials, each one contributes at its own power of `g`. A nonzero `ε_l` needs `l` to be a sum of active indices `p − 2`, repeats allowed. In particular, the gcd of the active indices must divide `l`. That condition is necessary but not sufficient, and the first nonzero index alone does not decide which orders vanish.

For V = x²/2 + x⁶ + x⁸ the indices are 4 and 6. `l = 6` is not a multiple of 4, yet the x⁸ term contributes there at first order. And gcd(4, 6) = 2 divides `l = 2`, but `ε_2` vanishes because 2 is not a sum of 4s and 6s.

```julia
potm = Potential([1//2, 0//1, 0//1, 0//1, 1//1, 0//1, 1//1])
ε_l(potm, 0, 2)   # → 0//1
ε_l(potm, 0, 4)   # → 15//8    (⟨0|x⁶|0⟩, first order in x⁶)
ε_l(potm, 0, 6)   # → 105//16  (⟨0|x⁸|0⟩, first order in x⁸)
```

### Energy polynomial in ν

At each fixed perturbation order, the correction is a polynomial in ν. `find_epoly` fits that polynomial; `evaluate_epoly` evaluates it.

```julia
epoly = find_epoly(0, pot)        # → [0.5, 1.0]         (ε⁽⁰⁾(ν) = 0.5 + ν)
epoly = find_epoly(2, pot)        # → [0.75, 1.5, 1.5]   (ε⁽²⁾(ν) = 0.75 + 1.5ν + 1.5ν²)

evaluate_epoly(0, epoly)   # → 0.75   (matches ε_l(pot, 0, 2))
evaluate_epoly(1, epoly)   # → 3.75   (matches ε_l(pot, 1, 2))
evaluate_epoly(3, epoly)   # → 21.75
```

**Compute all orders up to 50:**

```julia
ε_polys = [find_epoly(n, pot) for n = 0:50]
```

**Derivatives of ε(ν) evaluated at ν = 0** — entry `k` is the k-th derivative:

```julia
ds = epoly_taylor_derivatives(find_epoly(2, pot))   # → [1.5, 3.0]
```

### Numeric precision modes

Three precision modes are supported; the output type matches the element type of `vcoeffs`.

**Float64** (default):

```julia
pot = Potential([0.5, 0.0, 1.0])
ε_l(pot, 0, 2)   # → 0.75
```

**Exact rational arithmetic** — pass `Rational` coefficients. Integer-typed rationals are automatically promoted to `Rational{BigInt}` to prevent overflow at high orders. The leading coefficient must satisfy 2·vcoeffs[1] being a perfect square so that ω is rational.

```julia
pot_r = Potential([1//2, 0//1, 1//1])
ε_l(pot_r, 0, 2)   # → 3//4   (exact)
ε_l(pot_r, 1, 2)   # → 15//4  (exact)

find_epoly(4, pot_r)   # → [-21//8, -59//8, -51//8, -17//4]
```

**Arbitrary precision** — pass `BigFloat` coefficients:

```julia
pot_bf = Potential(BigFloat.([0.5, 0.0, 1.0]))
ε_polys_bf = [find_epoly(n, pot_bf) for n = 0:50]
```

Float64, BigFloat, and Rational potentials each carry independent caches; no manual flushing is needed.

### Iterative (type-stable) API

For computing all orders 0…maxorder at a fixed ν, `fill_Akl!` avoids the overhead of the recursive cache lookups:

```julia
ν, maxorder = 2, 10
Akl, ε = initialize_Akl_eps(pot, ν, maxorder)
fill_Akl!(Akl, ε, pot, ν, maxorder)
# ε[l+1] now equals ε_l(pot, ν, l) for every l in 0:maxorder
```

**Which API to use?** Reach for `fill_Akl!` whenever you need many orders at a fixed ν — it is type-stable and substantially faster than repeated `ε_l` calls. Use the recursive `ε_l` / `A_kl` for ad-hoc single-value queries or when you do not know up front how many orders you will need; results are memoized inside `pot` and reused across calls. For very high orders (`l` in the hundreds, plausible in BigFloat asymptotic studies) prefer `fill_Akl!` — the recursive path can exhaust Julia's default stack at that depth.

### Wave-function coefficients

If you only need the eigenstate expansion (no energy-correction array), `eigenstate_coeffs` is a thin wrapper around `initialize_Akl_eps` + `fill_Akl!` that returns just the 2D coefficient matrix:

```julia
Akl = eigenstate_coeffs(pot, 0, 4)   # ground state coefficients up to BW order 4
Akl[1, 1] == 1.0                      # k = ν normalisation at order 0
```

The matrix is indexed as `Akl[k+1, l+1]` (1-based offset), where `k` is the BW expansion index and `l` is the BW recursion order. The ν-th perturbed eigenstate is ψ_ν(x) = e^{-ω x² / 2}·∑_l g^l ∑_k Akl[k+1, l+1] · x^k, with g the coupling of H(g) and ω = √(2·vcoeffs[1]).

## API reference

| Function | Description |
|---|---|
| `Potential(vcoeffs)` | Construct a potential from a coefficient vector; owns memoization caches |
| `clear_cache!(pot)` | Empty the memoization caches inside `pot`; useful for long sweeps |
| `ε_l(pot, ν, l)` | Perturbative energy correction at order `l` for quantum number `ν` |
| `A_kl(pot, ν, k, l)` | Wave function expansion coefficient A_{k,l}^(ν) |
| `max_k(pot, ν, l)` | Upper bound on the k-index at perturbation order `l` |
| `find_epoly(order, pot)` | Fit the order-`order` energy correction as a polynomial in ν; returns coefficient vector |
| `evaluate_epoly(n, epoly)` | Evaluate an energy polynomial at ν = n |
| `epoly_taylor_derivatives(epoly)` | Derivatives of an energy polynomial evaluated at ν = 0 (entry `k` is the k-th derivative) |
| `initialize_Akl_eps(pot, ν, l)` | Allocate zero arrays for the iterative solver |
| `fill_Akl!(Akl, ε, pot, ν, maxorder)` | Fill pre-allocated arrays in-place (iterative, type-stable) |
| `eigenstate_coeffs(pot, ν, l)` | Wave function expansion coefficients up to BW recursion order `l` (2D array) |

## Tests

```julia
using Pkg; Pkg.test()
```

## Contributing

Use one branch per feature or fix and open a pull request into `main`.
See [CONTRIBUTING.md](CONTRIBUTING.md) for CI checks and the release process.

## Benchmarks

End-to-end timings against the reference Mathematica implementation `BenderWu.m` for both Float64 and exact-rational arithmetic are reported in [BENCHMARKS.md](BENCHMARKS.md).
See [benchmark/README.md](benchmark/README.md) for how to reproduce them.
