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

Dependencies: `Memoize` (runtime), `BenchmarkTools` (dev/notebook only).

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

All code lives in `BenderWu.jl`. The algorithm computes perturbative corrections to energy eigenvalues $E_\nu$ of a Hamiltonian with potential given by `vcoeffs` (a vector where `vcoeffs[n]` is the coefficient of $x^n$ in the potential).

### Key Abstractions

- **`vcoeffs`** — potential polynomial coefficients, e.g. `[0.5, 0.0, 1.0]` for $V = 0.5x^0 + x^2$. The first coefficient sets $\omega = \sqrt{2 \cdot \text{vcoeffs}[1]}$.
- **`ν`** — quantum number (energy level index, 0-based)
- **`l`** — perturbation order (only even orders contribute; odd orders return zero)

### Recursive (memoized) implementation

- `max_k(ν, l, vcoeffs)` — upper bound on k-index at perturbation order `l`
- `A_kl(ν, k, l, vcoeffs)` — wave function expansion coefficient; memoized with `@memoize`, mutually recursive with `ε_l`
- `ε_l(ν, l, vcoeffs)` — energy correction at order `l`; memoized, mutually recursive with `A_kl`

`@memoize` caches calls globally per Julia session. Call `Memoize.empty_all_caches!()` when changing `vcoeffs` or switching precision (Float64 vs BigFloat).

### Iterative (type-stable) implementation

- `initialize_Akl_eps(ν, l, vcoeffs)` — pre-allocates arrays; output type is inferred from `vcoeffs[1]`
- `fill_Akl!(Akl, ε, ν, maxorder, vcoeffs)` — fills arrays in-place following three steps per order: compute `Akl` for k > ν, compute `ε`, compute `Akl` for k < ν

### Energy polynomial fitting

- `find_epoly(order, vcoeffs)` — solves a linear system to express $\varepsilon(\nu)$ at perturbation `order` as a polynomial in $\nu$; returns coefficient vector
- `find_epoly_derivative(epoly)` — differentiates the polynomial (uses `big` factorial for orders ≥ 20)
- `evaluate_epoly(n, epoly)` — evaluates the polynomial at a given $\nu = n$

### Precision

Pass `BigFloat` coefficients (e.g. `[big(0.5), big(0.0), big(1.0)]`) to get arbitrary-precision results throughout all functions. The output type is always inferred from `vcoeffs[1]`.
