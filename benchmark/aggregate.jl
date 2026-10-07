#!/usr/bin/env julia
# Loads benchmark/results_julia.json and benchmark/results_mma.json, validates
# that both sides ran the same cases and agree exactly on ε_l, draws plots, and
# writes BENCHMARKS.md at the repo root. Any validation failure throws before
# anything is written.
#
# Usage:
#   julia --project=benchmark benchmark/aggregate.jl

using JSON3
using Plots
using Printf

include(joinpath(@__DIR__, "cases.jl"))

const ROOT = abspath(joinpath(@__DIR__, ".."))
const PLOTS_DIR = joinpath(@__DIR__, "plots")
isdir(PLOTS_DIR) || mkpath(PLOTS_DIR)

# ---------------------------------------------------------------------------
# Loading

load_results(path) = JSON3.read(read(path, String))

# ---------------------------------------------------------------------------
# Validation (exact rationals — bit-for-bit equality)

"""
Parse a serialised exact ε value (Julia "a//b" or Mathematica "a/b", possibly
parenthesised) into a Rational{BigInt}.
"""
function parse_eps(s::AbstractString)
    s = strip(strip(s), ['(', ')'])
    if occursin("//", s)
        num, den = split(s, "//")
        return Rational{BigInt}(parse(BigInt, strip(num)),
                                parse(BigInt, strip(den)))
    elseif occursin("/", s)
        num, den = split(s, "/")
        return Rational{BigInt}(parse(BigInt, strip(num)),
                                parse(BigInt, strip(den)))
    else
        return Rational{BigInt}(parse(BigInt, s))
    end
end

"""
Check one cell's ε vectors and return a list of problems (empty if they agree).
N is the highest power of g. Julia stores ε_0, ε_1, …, ε_N with odd-l entries
zero; Mathematica stores ε_0, ε_2, …, ε_N.
"""
function eps_problems(jl_eps::Vector, mma_eps::Vector, N::Int)
    problems = String[]
    length(jl_eps) == N + 1 ||
        push!(problems, "Julia has $(length(jl_eps)) entries, expected $(N + 1)")
    length(mma_eps) == N ÷ 2 + 1 ||
        push!(problems, "Mathematica has $(length(mma_eps)) entries, expected $(N ÷ 2 + 1)")
    isempty(problems) || return problems
    all(iszero, jl_eps[2:2:end]) ||
        push!(problems, "Julia has nonzero odd-order entries")
    jl_eps[1:2:end] == mma_eps ||
        push!(problems, "even-order ε_l differ")
    return problems
end

# ---------------------------------------------------------------------------
# Joining

const Key = Tuple{String, Int, Int}

"""
Index one side's results by (potential, ν, N), throwing on duplicate keys.
"""
function index_results(results, side)
    by = Dict{Key, Any}()
    for r in results
        k = (String(r.potential), Int(r.nu), Int(r.N))
        haskey(by, k) && error("duplicate $side result for $k")
        by[k] = r
    end
    return by
end

"""
Join Julia and Mathematica results cell by cell. Throws if either side is
missing a case, if there are no cases, or if any ε vector disagrees.
"""
function join_results(jl, mma)
    jl_by  = index_results(jl.results, "Julia")
    mma_by = index_results(mma.results, "Mathematica")

    only_jl  = sort!(collect(setdiff(keys(jl_by), keys(mma_by))))
    only_mma = sort!(collect(setdiff(keys(mma_by), keys(jl_by))))
    if !isempty(only_jl) || !isempty(only_mma)
        error("result sets differ.\n  missing in Mathematica: $only_jl\n",
              "  missing in Julia: $only_mma")
    end
    isempty(jl_by) && error("no benchmark results to compare")

    rows = []
    failures = String[]
    for k in sort!(collect(keys(jl_by)))
        j, m = jl_by[k], mma_by[k]
        jl_eps  = [parse_eps(s) for s in j.epsilons]
        mma_eps = [parse_eps(s) for s in m.epsilons]
        for p in eps_problems(jl_eps, mma_eps, k[3])
            push!(failures, "$(k[1]) ν=$(k[2]) N=$(k[3]): $p")
        end
        push!(rows, (
            potential = k[1],
            nu = k[2],
            N = k[3],
            jl_ms  = Float64(j.time_ns_median) / 1e6,
            mma_ms = Float64(m.time_ns_median) / 1e6,
            ratio  = Float64(m.time_ns_median) / Float64(j.time_ns_median),
        ))
    end
    if !isempty(failures)
        foreach(f -> println("  MISMATCH: ", f), failures)
        error("$(length(failures)) validation failure(s); nothing written")
    end
    return rows
end

# ---------------------------------------------------------------------------
# Plotting

# Potential names in POTENTIALS order, restricted to those present in rows.
potential_names(rows) = [p.name for p in POTENTIALS if any(r -> r.potential == p.name, rows)]

function plot_per_potential(rows)
    paths = String[]
    for pot in potential_names(rows)
        sub = filter(r -> r.potential == pot, rows)
        nus = sort!(unique(r.nu for r in sub))
        plt = plot(; xlabel = "maximum power of g, N",
                   ylabel = "median time (ms)",
                   yscale = :log10,
                   title  = "$pot — exact rational arithmetic",
                   legend = :topleft,
                   size   = (700, 450))
        for ν in nus
            rs = sort(filter(r -> r.nu == ν, sub), by = r -> r.N)
            Ns = [r.N for r in rs]
            plot!(plt, Ns, [r.jl_ms  for r in rs];
                  marker = :circle, label = "Julia ν=$ν")
            plot!(plt, Ns, [r.mma_ms for r in rs];
                  marker = :square, linestyle = :dash, label = "Mathematica ν=$ν")
        end
        path = joinpath(PLOTS_DIR, "$(pot).png")
        savefig(plt, path)
        push!(paths, path)
    end
    return paths
end

function plot_speedup(rows)
    plt = plot(; xlabel = "maximum power of g, N",
               ylabel = "Mathematica / Julia (median time ratio)",
               yscale = :log10,
               title  = "Speedup factor — exact rational arithmetic",
               legend = :topright,
               size   = (800, 500))
    for pot in potential_names(rows)
        rs = sort(filter(r -> r.potential == pot && r.nu == 0, rows),
                  by = r -> r.N)
        isempty(rs) && continue
        Ns = [r.N for r in rs]
        plot!(plt, Ns, [r.ratio for r in rs];
              marker = :circle, label = pot)
    end
    hline!(plt, [1.0]; color = :gray, linestyle = :dot, label = "parity")
    path = joinpath(PLOTS_DIR, "speedup.png")
    savefig(plt, path)
    return path
end

# ---------------------------------------------------------------------------
# Markdown report

function fmt_ms(x)
    x < 0.01 && return @sprintf("%.4f", x)
    x < 1.0  && return @sprintf("%.3f", x)
    x < 100  && return @sprintf("%.2f", x)
    return @sprintf("%.1f", x)
end

function fmt_ratio(x)
    x < 1   && return @sprintf("%.2f×", x)
    x < 100 && return @sprintf("%.1f×", x)
    return @sprintf("%.0f×", x)
end

function write_report(rows, jl_machine, mma_machine; outpath)
    pots = potential_names(rows)

    io = IOBuffer()
    println(io, "# BenderWu — Julia vs Mathematica benchmarks\n")
    println(io, "Comparison of this Julia package against the reference Mathematica")
    println(io, "implementation `BenderWu.m` ([arXiv:1608.08256](https://arxiv.org/abs/1608.08256)). Both")
    println(io, "implementations compute the perturbative energy corrections ε_l of")
    println(io, "E = Σ ε_l g^l at fixed quantum number ν for the polynomial potentials")
    println(io, "shown below. N is the highest power of g: Julia computes l = 0…N, and")
    println(io, "Mathematica is called as `BenderWu[V, x, ν, N/2]`, since its order")
    println(io, "argument counts powers of g².\n")

    println(io, "Only **exact-rational arithmetic** is benchmarked.")
    println(io, "Mathematica's `BenderWu` evaluates the recursion through its")
    println(io, "symbolic term-rewriting pipeline regardless of coefficient")
    println(io, "precision, so the gap there mostly measures evaluator overhead")
    println(io, "rather than algorithmic efficiency. In exact-rational mode both")
    println(io, "sides do genuine big-integer arithmetic and the comparison is")
    println(io, "meaningful.\n")

    println(io, "## Setup\n")
    println(io, "| | |")
    println(io, "|---|---|")
    println(io, "| Julia      | $(jl_machine["julia_version"]) |")
    println(io, "| Mathematica| $(mma_machine["mathematica_version"]) |")
    println(io, "| CPU        | $(jl_machine["cpu"]) ($(jl_machine["ncores"]) threads) |")
    println(io, "| OS         | $(jl_machine["os"]) / $(jl_machine["arch"]) |")
    println(io)

    println(io, "Both sides report the median over up to 50 samples within a 5 s")
    println(io, "budget per case, after a warm-up call. Julia uses")
    println(io, "`BenchmarkTools.@benchmark` with a fresh `Potential` per sample so")
    println(io, "caches are cold; Mathematica uses `AbsoluteTiming` in a loop with")
    println(io, "the same limits.\n")

    println(io, "## Validation\n")
    println(io, "For every case, Julia's ε vector has N + 1 entries with all odd")
    println(io, "orders zero, Mathematica's has N/2 + 1 entries, and after dropping")
    println(io, "Julia's odd orders the two vectors are identical — bit-for-bit")
    println(io, "equality on rationals. Both sides ran the same set of cases.\n")
    println(io, "**$(length(rows))/$(length(rows))** cases match exactly.\n")

    println(io, "## Potentials\n")
    println(io, "| name | V(x) |")
    println(io, "|---|---|")
    for p in POTENTIALS
        p.name in pots && println(io, "| $(p.name) | $(p.display) |")
    end
    println(io)

    println(io, "## Results\n")
    println(io, "### Per-potential timings\n")
    for pot in pots
        img = "benchmark/plots/$(pot).png"
        isfile(joinpath(ROOT, img)) || continue
        println(io, "![$pot]($img)\n")
    end

    speedup_path = "benchmark/plots/speedup.png"
    if isfile(joinpath(ROOT, speedup_path))
        println(io, "### Speedup factor\n")
        println(io, "Per-cell ratio of Mathematica median time to Julia median time")
        println(io, "(at ν = 0). Higher is better for Julia.\n")
        println(io, "![Speedup]($speedup_path)\n")
    end

    println(io, "### Timing table\n")
    println(io, "| potential | ν | N | Julia | Mathematica | speedup |")
    println(io, "|---|---|---|---:|---:|---:|")
    order = Dict(name => i for (i, name) in enumerate(pots))
    for r in sort(rows, by = r -> (order[r.potential], r.nu, r.N))
        println(io, "| $(r.potential) | $(r.nu) | $(r.N) | ",
                fmt_ms(r.jl_ms), " ms | ",
                fmt_ms(r.mma_ms), " ms | ",
                fmt_ratio(r.ratio), " |")
    end
    println(io)

    println(io, "## Reproducing\n")
    println(io, "See [benchmark/README.md](benchmark/README.md) for the exact")
    println(io, "commands to regenerate this report.")

    open(outpath, "w") do f
        write(f, take!(io))
    end
end

# ---------------------------------------------------------------------------
# Main

function main()
    jl  = load_results(joinpath(@__DIR__, "results_julia.json"))
    mma = load_results(joinpath(@__DIR__, "results_mma.json"))

    rows = join_results(jl, mma)
    println("All $(length(rows)) cases match exactly.")

    plot_per_potential(rows)
    plot_speedup(rows)

    jl_machine  = Dict(string(k) => v for (k, v) in pairs(jl.machine))
    mma_machine = Dict(string(k) => v for (k, v) in pairs(mma.machine))

    outpath = joinpath(ROOT, "BENCHMARKS.md")
    write_report(rows, jl_machine, mma_machine; outpath)
    println("Wrote ", outpath)
end

main()
