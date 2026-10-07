# AGENTS.md

## Project Overview

Julia implementation of the **Bender-Wu method** for computing perturbative energy levels of quantum systems with polynomial potentials (reference: arXiv:1608.08256). The project is a single-file implementation (`BenderWu.jl`).

## Environment Setup

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using BenderWu
```

Dependencies: none.

## Running and Testing

Run the test suite:

```julia
using Pkg; Pkg.test()
```

Or directly:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Architecture

All code lives in `src/BenderWu.jl`. The algorithm computes perturbative corrections to energy eigenvalues $E_\nu$ of a Hamiltonian with potential given by `vcoeffs` (a vector where `vcoeffs[n]` is the coefficient of $x^{n+1}$ in the potential).

### Key Abstractions

- **`Potential{T}`** — wraps `vcoeffs` (potential polynomial coefficients) and owns its own memoization caches. Create one per potential and reuse it. Cache lifetime is GC-managed; call `clear_cache!(pot)` to drop entries early during long sweeps. Cache access is guarded by a `ReentrantLock`, so a `Potential` can be shared across threads.
  - `vcoeffs[n]` is the coefficient of $x^{n+1}$, e.g. `Potential([0.5, 0.0, 1.0])` for $V = \frac{1}{2}x^2 + x^4$
  - The first coefficient sets $\omega = \sqrt{2 \cdot \text{vcoeffs}[1]}$
- **`ν`** — quantum number (energy level index, 0-based)
- **`l`** — perturbation order (only even orders contribute; odd orders return zero)

### Precision and numeric types

The output type is always inferred from `eltype(pot.vcoeffs)`. Three modes are supported:

- **Float64** (default): `Potential([0.5, 0.0, 1.0])`
- **BigFloat** (arbitrary precision): `Potential(BigFloat.([0.5, 0.0, 1.0]))`
- **Rational** (exact arithmetic): `Potential([1//2, 0//1, 1//1])` — integer-typed rationals (`Rational{Int64}` etc.) are automatically promoted to `Rational{BigInt}` to prevent overflow at higher perturbation orders.

For rational potentials, `_compute_ω` is dispatched to an exact method that uses `isqrt` and validates that $2 \cdot \text{vcoeffs}[1]$ is a perfect square. Float64 and BigFloat `Potential` objects have fully independent caches; no manual flushing is needed.

## Conventions

### GitHub issues

* Check open and closed issues for duplicates before filing; link related issues.
* Keep each issue focused on one independently actionable outcome. Track broader design proposals separately from fixes.
* Use concise, conventional titles such as `fix: make solver buffers safe to reuse` or `perf: skip coefficients that vanish by symmetry`.
* Describe the concrete problem and its impact, then include a minimal runnable reproduction with actual and expected results. Record the affected commit, Julia version, and relevant numeric types and inputs.
* For performance findings, include measured timings, equivalent workloads, hardware, warm-up, and sampling method. Identify prototype results and indicative measurements as such, and include enough code to reproduce them.
* Define a clear, verifiable expected outcome. Separate confirmed behavior from proposed solutions and unresolved design choices.
* Link to relevant source using commit permalinks so the evidence remains readable after the code changes.

### Pull requests and releases

* Start each feature or fix on its own branch from the latest `main`.
* Open draft pull requests targeting `main`; do not use `dev` as an integration branch.
* Rebase onto the latest `main` before opening a pull request when explicitly authorized to rebase.
* Use conventional commit style for pull request titles.
* Describe the problem, solution, and validation concisely.
* Disclose substantial agent-generated changes in the pull request description.
* Follow `CONTRIBUTING.md` for CI and release preparation.
* Release version changes and `.github/release-notes/vX.Y.Z.md` go through a pull request into `main` before tagging.

### Git commits

Use conventional commits: `<type>[optional scope]: <description>`

Examples:

* `feat: add CSV export`
* `fix(api): handle empty responses`
* `docs: update setup guide`
* `refactor(parser): simplify token handling`
* `test: cover invalid input`
* `chore: update dependencies`

Allowed types: `feat`, `fix`, `docs`, `refactor`, `test`, `perf`, `build`, `ci`, `style`, `chore`, `revert`.

Rules:

* Use imperative, lowercase descriptions.
* Do not end the subject with a period.
* Keep the subject concise, preferably under 72 characters.
* Make each commit one logical change.
* Use `!` for breaking changes: `feat(api)!: replace pagination`.
* Add a body when the reason or impact is not obvious.
* Do not commit, amend, rebase, squash, or force-push unless explicitly requested.
* Do not commit in plans. Plans stop at written-and-verified.
* Do not add `Co-authored-by` or other AI-attribution trailers unless explicitly required by the repository or requested by the user.
