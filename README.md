# BenderWu.jl

Julia implementation of the Bender-Wu method for computing perturbative energy levels of quantum systems with polynomial potentials, following the recursive relations in [arXiv:1608.08256](https://arxiv.org/abs/1608.08256).

## Setup

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using BenderWu
```

## Usage

A potential is represented as a coefficient vector `vcoeffs`, where `vcoeffs[n]` is the coefficient of $x^{n-1}$. For example, `[0.5, 0.0, 1.0]` encodes $V(x) = \frac{1}{2}x^2 + x^4$ (quartic oscillator with $\omega = 1$).

**Energy correction at a given perturbation order:**

```julia
vcoeffs = [0.5, 0.0, 1.0]

ε_l(0, 0, vcoeffs)   # zeroth-order energy of ground state: 0.5
ε_l(1, 2, vcoeffs)   # second-order correction for ν=1: 3.75
```

Only even orders contribute; odd orders return zero exactly.

**Fit the energy as a polynomial in the quantum number ν:**

```julia
epoly = find_epoly(2, vcoeffs)       # coefficients of ε(ν) at order l=2
evaluate_epoly(3, epoly)             # evaluate at ν=3
```

**Compute all orders up to 50:**

```julia
ε_polys = [find_epoly(n, vcoeffs) for n = 0:50]
```

**Arbitrary precision** — pass `BigFloat` coefficients:

```julia
vcoeffs_bf = BigFloat.([0.5, 0.0, 1.0])
ε_polys_bf = [find_epoly(n, vcoeffs_bf) for n = 0:50]
```

**Derivative of the energy polynomial** (Taylor coefficients at ν=0):

```julia
find_epoly_derivative(epoly)
```

**Iterative (non-memoized) alternative:**

```julia
ν, maxorder = 2, 10
Akl, ε = initialize_Akl_eps(ν, maxorder, vcoeffs)
fill_Akl!(Akl, ε, ν, maxorder, vcoeffs)
```

### Memoization cache

The memoized functions (`ε_l`, `A_kl`, `max_k`) cache results globally. Clear the cache when switching between different `vcoeffs` inputs or between Float64 and BigFloat to avoid stale entries:

```julia
using Memoize
Memoize.empty_all_caches!()
```

## Tests

```julia
using Pkg; Pkg.test()
```
