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

A potential is represented by a `Potential` object constructed from a coefficient vector, where `vcoeffs[n]` is the coefficient of $x^{n-1}$. For example, `[0.5, 0.0, 1.0]` encodes $V(x) = \frac{1}{2}x^2 + x^4$ (quartic oscillator with $\omega = 1$).

```julia
pot = Potential([0.5, 0.0, 1.0])
```

Create one `Potential` per potential and reuse it — results are cached inside the struct and freed automatically when it goes out of scope.

**Energy correction at a given perturbation order:**

```julia
ε_l(pot, 0, 0)   # zeroth-order energy of ground state: 0.5
ε_l(pot, 1, 2)   # second-order correction for ν=1: 3.75
```

Only even orders contribute; odd orders return zero exactly.

**Fit the energy as a polynomial in the quantum number ν:**

```julia
epoly = find_epoly(2, pot)       # coefficients of ε(ν) at order l=2
evaluate_epoly(3, epoly)         # evaluate at ν=3
```

**Compute all orders up to 50:**

```julia
ε_polys = [find_epoly(n, pot) for n = 0:50]
```

**Arbitrary precision** — pass `BigFloat` coefficients:

```julia
pot_bf = Potential(BigFloat.([0.5, 0.0, 1.0]))
ε_polys_bf = [find_epoly(n, pot_bf) for n = 0:50]
```

Float64 and BigFloat potentials have independent caches; no manual flushing needed.

**Derivative of the energy polynomial** (Taylor coefficients at ν=0):

```julia
find_epoly_derivative(epoly)
```

**Iterative (non-recursive) alternative:**

```julia
ν, maxorder = 2, 10
Akl, ε = initialize_Akl_eps(pot, ν, maxorder)
fill_Akl!(Akl, ε, pot, ν, maxorder)
```

## Tests

```julia
using Pkg; Pkg.test()
```
