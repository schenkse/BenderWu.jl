# BenderWu.jl

Julia implementation of the **Bender-Wu method** for computing perturbative energy corrections to eigenvalues of 1D quantum systems with polynomial potentials.

The Hamiltonian is H = p²/2 + V(x), where V(x) = Σ vcoeffs[n] · xⁿ⁺¹. The energy eigenvalue E_ν is expanded order by order in a coupling constant; this package computes those perturbative corrections using the Bender-Wu recursive relations.

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

The frequency ω is derived automatically from the leading term: ω = √(2·vcoeffs[1]).

Create one `Potential` per potential and reuse it — results are memoized inside the struct and freed automatically when it goes out of scope.

### Energy corrections ε_l(pot, ν, l)

`ε_l(pot, ν, l)` returns the perturbative energy correction at order `l` for quantum number `ν`.

```julia
ε_l(pot, 0, 0)   # → 0.5    (ground state, unperturbed: ω·(0 + 1/2))
ε_l(pot, 1, 0)   # → 1.5    (ν=1, unperturbed: ω·(1 + 1/2))
ε_l(pot, 0, 2)   # → 0.75   (first non-trivial correction for ν=0)
ε_l(pot, 1, 2)   # → 3.75
ε_l(pot, 2, 2)   # → 9.75
```

Odd-order corrections vanish identically and are returned as exact zero.

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

## API reference

| Function | Description |
|---|---|
| `Potential(vcoeffs)` | Construct a potential from a coefficient vector; owns memoization caches |
| `ε_l(pot, ν, l)` | Perturbative energy correction at order `l` for quantum number `ν` |
| `A_kl(pot, ν, k, l)` | Wave function expansion coefficient A_{k,l}^(ν) |
| `max_k(pot, ν, l)` | Upper bound on the k-index at perturbation order `l` |
| `find_epoly(order, pot)` | Fit the order-`order` energy correction as a polynomial in ν; returns coefficient vector |
| `evaluate_epoly(n, epoly)` | Evaluate an energy polynomial at ν = n |
| `epoly_taylor_derivatives(epoly)` | Derivatives of an energy polynomial evaluated at ν = 0 (entry `k` is the k-th derivative) |
| `initialize_Akl_eps(pot, ν, l)` | Allocate zero arrays for the iterative solver |
| `fill_Akl!(Akl, ε, pot, ν, maxorder)` | Fill pre-allocated arrays in-place (iterative, type-stable) |

## Tests

```julia
using Pkg; Pkg.test()
```
