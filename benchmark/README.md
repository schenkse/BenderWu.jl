# Benchmarks

Compares this Julia BenderWu implementation against the reference Mathematica
package `BenderWu.m` by Sulejmanpasic (the implementation released with
arXiv:1608.08256). The headline numbers and plots live in
[../BENCHMARKS.md](../BENCHMARKS.md).

## Files

| File | Purpose |
|---|---|
| [`cases.jl`](cases.jl) | Test matrix — potentials, ν values, order sweeps. Single source of truth shared by the Julia driver and the aggregator. |
| [`run_julia.jl`](run_julia.jl) | Times Julia via `BenchmarkTools.@benchmark`. Writes `results_julia.json`. |
| [`run_mathematica.wls`](run_mathematica.wls) | Times Mathematica via `RepeatedTiming`. Writes `results_mma.json`. |
| [`aggregate.jl`](aggregate.jl) | Validates ε_l agreement, draws plots into `plots/`, and writes `../BENCHMARKS.md`. |
| `Project.toml` | Julia environment for the benchmark scripts (kept separate so the BenderWu package itself stays dependency-free). |

## What is measured

Both implementations compute the perturbative energy corrections ε_l for
l = 0…N at fixed quantum number ν, given a polynomial potential V(x). The
Julia path uses `initialize_Akl_eps` + `fill_Akl!`. The Mathematica path
uses `BenderWu[V, x, ν, N, Output -> "Energy", OutputStyle -> "Array"]`.

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
- Mathematica with `wolframscript` (the script falls back to the macOS
  default location `/Applications/Wolfram.app/Contents/MacOS/wolframscript`
  if `wolframscript` is not on `PATH`)
- The Mathematica package `BenderWu.m` installed at
  `~/Library/Wolfram/Applications/BenderWu.m` (Sulejmanpasic, available
  alongside arXiv:1608.08256)

One-time setup of the benchmark environment:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

Run the full sweeps and regenerate `BENCHMARKS.md`:

```bash
# Julia (~30 s)
julia --project=benchmark benchmark/run_julia.jl

# Mathematica (~10 min)
wolframscript -file benchmark/run_mathematica.wls
# or, if wolframscript is not on PATH:
/Applications/Wolfram.app/Contents/MacOS/wolframscript -file benchmark/run_mathematica.wls

# Aggregate (validate + plot + write report)
julia --project=benchmark benchmark/aggregate.jl
```

For a fast sanity check, both runners support `--quick` which restricts the
sweep to a single ν and a single order per mode:

```bash
julia --project=benchmark benchmark/run_julia.jl --quick
wolframscript -file benchmark/run_mathematica.wls --quick
julia --project=benchmark benchmark/aggregate.jl
```

## Methodology notes

- Each Julia sample rebuilds the `Potential` (and thus its cache) so caches
  are cold. This matches Mathematica, where every `BenderWu[…]` call
  recomputes from scratch.
- Mathematica's `RepeatedTiming` chooses the number of repetitions
  automatically — many for fast cells, a single run for slow ones — and
  returns the average.
- Validation requires bit-for-bit equality of the rational ε_l vectors.
- Mathematica's per-call setup cost (Module initialisation, OptionsPattern
  parsing, symbolic preprocessing) is included in its timings — these are
  end-to-end "compute one sweep" numbers, not isolated kernel timings.
