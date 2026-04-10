# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Julia implementation of the **Bender-Wu method** for computing perturbative energy levels of quantum systems with polynomial potentials (reference: arXiv:1608.08256). The project is a single-file implementation (`BenderWu.jl`) with an interactive Jupyter notebook for exploration.

## Environment Setup

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using BenderWu
```

Dependencies: none at runtime. `BenchmarkTools` is used in the notebook only.

## Running and Testing

Run the test suite:

```julia
using Pkg; Pkg.test()
```

Or directly:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Interactive exploration via Jupyter:

```bash
jupyter notebook BenderWu_Interactive.ipynb
```

## Architecture

All code lives in `src/BenderWu.jl`. The algorithm computes perturbative corrections to energy eigenvalues $E_\nu$ of a Hamiltonian with potential given by `vcoeffs` (a vector where `vcoeffs[n]` is the coefficient of $x^{n-1}$ in the potential).

### Key Abstractions

- **`Potential{T}`** — wraps `vcoeffs` (potential polynomial coefficients) and owns its own memoization caches. Create one per potential and reuse it. Cache lifetime is GC-managed.
  - `vcoeffs[n]` is the coefficient of $x^{n-1}$, e.g. `Potential([0.5, 0.0, 1.0])` for $V = \frac{1}{2}x^2 + x^4$
  - The first coefficient sets $\omega = \sqrt{2 \cdot \text{vcoeffs}[1]}$
- **`ν`** — quantum number (energy level index, 0-based)
- **`l`** — perturbation order (only even orders contribute; odd orders return zero)

### Recursive (cached) implementation

- `max_k(pot, ν, l)` — upper bound on k-index at perturbation order `l`
- `A_kl(pot, ν, k, l)` — wave function expansion coefficient; cached in `pot._Akl_cache`, mutually recursive with `ε_l`
- `ε_l(pot, ν, l)` — energy correction at order `l`; cached in `pot._εl_cache`, mutually recursive with `A_kl`

### Iterative (type-stable) implementation

- `initialize_Akl_eps(pot, ν, l)` — pre-allocates arrays; output type is inferred from `eltype(pot.vcoeffs)`
- `fill_Akl!(Akl, ε, pot, ν, maxorder)` — fills arrays in-place following three steps per order: compute `Akl` for k > ν, compute `ε`, compute `Akl` for k < ν

### Energy polynomial fitting

- `find_epoly(order, pot)` — solves a linear system to express $\varepsilon(\nu)$ at perturbation `order` as a polynomial in $\nu$; returns coefficient vector
- `find_epoly_derivative(epoly)` — differentiates the polynomial (uses `big` factorial for orders ≥ 20)
- `evaluate_epoly(n, epoly)` — evaluates the polynomial at a given $\nu = n$

### Precision and numeric types

The output type is always inferred from `eltype(pot.vcoeffs)`. Three modes are supported:

- **Float64** (default): `Potential([0.5, 0.0, 1.0])`
- **BigFloat** (arbitrary precision): `Potential(BigFloat.([0.5, 0.0, 1.0]))`
- **Rational** (exact arithmetic): `Potential([1//2, 0//1, 1//1])` — integer-typed rationals (`Rational{Int64}` etc.) are automatically promoted to `Rational{BigInt}` to prevent overflow at higher perturbation orders.

For rational potentials, `_compute_ω` is dispatched to an exact method that uses `isqrt` and validates that $2 \cdot \text{vcoeffs}[1]$ is a perfect square. Float64 and BigFloat `Potential` objects have fully independent caches; no manual flushing is needed.
