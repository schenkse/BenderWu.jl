# Benchmarks

Compares this Julia BenderWu implementation against the reference Mathematica package `BenderWu.m` (the implementation released with [arXiv:1608.08256](https://arxiv.org/abs/1608.08256)). The headline numbers and plots live in [../BENCHMARKS.md](../BENCHMARKS.md).

## Files

| File | Purpose |
|---|---|
| [`cases.jl`](cases.jl) | Test matrix — potentials, ν values, order sweeps. Single source of truth shared by the Julia driver and the aggregator. |
| [`run_julia.jl`](run_julia.jl) | Times Julia via `BenchmarkTools.@benchmark`. Writes `results_julia.json`. |
| [`run_mathematica.wls`](run_mathematica.wls) | Times Mathematica via an `AbsoluteTiming` loop with the same limits as the Julia runner. Writes `results_mma.json`. |
| [`aggregate.jl`](aggregate.jl) | Validates ε_l agreement, draws plots into `plots/`, and writes `../BENCHMARKS.md`. |
| `Project.toml` | Julia environment for the benchmark scripts (kept separate so the BenderWu package itself stays dependency-free). |

## What is measured

Both implementations compute the perturbative energy corrections ε_l of
E = Σ ε_l g^l at fixed quantum number ν, given a polynomial potential V(x).
N is the highest power of g and must be even. The Julia path uses
`initialize_Akl_eps` + `fill_Akl!` for l = 0…N. The Mathematica path uses
`BenderWu[V, x, ν, N/2, Output -> "Energy", OutputStyle -> "Array"]`,
because its order argument counts powers of g².

Only **exact-rational arithmetic** is benchmarked: Julia
`Rational{BigInt}` against Mathematica's exact symbolic arithmetic. A
Float64 vs MachinePrecision comparison would not be meaningful — the
Mathematica package evaluates its recursion through the symbolic
term-rewriting pipeline regardless of coefficient precision, so any
gap there reflects evaluator overhead rather than algorithmic
efficiency. In exact mode both sides do genuine big-integer arithmetic
and the comparison is apples-to-apples.

## Running

You need:

- Julia ≥ 1.10
- Mathematica with `wolframscript`
- The Mathematica package `BenderWu.m` installed

One-time setup of the benchmark environment:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

Run the full sweeps and regenerate `BENCHMARKS.md`:

```bash
# Julia
julia --project=benchmark benchmark/run_julia.jl

# Mathematica
wolframscript -file benchmark/run_mathematica.wls

# Aggregate (validate + plot + write report)
julia --project=benchmark benchmark/aggregate.jl
```

For a fast sanity check, both runners support `--quick` which restricts the
sweep to ν = 0 and N = 10 for each potential:

```bash
julia --project=benchmark benchmark/run_julia.jl --quick
wolframscript -file benchmark/run_mathematica.wls --quick
julia --project=benchmark benchmark/aggregate.jl
```

## Methodology notes

- Each Julia sample rebuilds the `Potential` (and thus its cache) so caches
  are cold. This matches Mathematica, where every `BenderWu[…]` call
  recomputes from scratch.
- Both runners take up to 50 samples within a 5 s budget per case (at
  least one), after a warm-up call, and report the median and minimum.
- `aggregate.jl` requires both sides to contain exactly the same cases.
  For each case, Julia's ε vector must have N + 1 entries with all odd
  orders zero, and Mathematica's must have N/2 + 1 entries. After
  dropping Julia's odd orders, the two vectors must be bit-for-bit equal.
  Any failure throws before plots or `BENCHMARKS.md` are written.
- Mathematica's per-call setup cost (Module initialisation, OptionsPattern
  parsing, symbolic preprocessing) is included in its timings — these are
  end-to-end "compute one sweep" numbers, not isolated kernel timings.
