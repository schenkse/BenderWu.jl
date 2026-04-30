#!/usr/bin/env julia
# Loads benchmark/results_julia.json and benchmark/results_mma.json, validates
# that ε_l agree exactly across implementations, draws plots, and writes
# BENCHMARKS.md at the repo root.
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
Compare two exact ε vectors on the union of even-l indices. The Mathematica
list is {ε_0, ε_2, ε_4, ...}; the Julia list is {ε_0, ε_1, ε_2, ...} with
odd-l entries zero. We compare jl_eps[2k+1] with mma_eps[k+1].
"""
function compare_eps(jl_eps::Vector, mma_eps::Vector)
    n = min(length(mma_eps), (length(jl_eps) - 1) ÷ 2 + 1)
    for k in 0:(n-1)
        if jl_eps[2k + 1] != mma_eps[k + 1]
            return false
        end
    end
    return true
end

# ---------------------------------------------------------------------------
# Joining

function join_results(jl, mma)
    key(r) = (String(r.potential), Int(r.nu), Int(r.N))
    mma_by = Dict(key(r) => r for r in mma.results)
    rows = []
    for j in jl.results
        m = get(mma_by, key(j), nothing)
        m === nothing && continue
        jl_eps  = [parse_eps(s) for s in j.epsilons]
        mma_eps = [parse_eps(s) for s in m.epsilons]
        matched = compare_eps(jl_eps, mma_eps)
        push!(rows, (
            potential = String(j.potential),
            nu = Int(j.nu),
            N = Int(j.N),
            jl_ms  = Float64(j.time_ns_median) / 1e6,
            mma_ms = Float64(m.time_ns_median) / 1e6,
            ratio  = Float64(m.time_ns_median) / Float64(j.time_ns_median),
            matched = matched,
        ))
    end
    return rows
end

# ---------------------------------------------------------------------------
# Plotting

function plot_per_potential(rows)
    paths = String[]
    for pot in unique(r.potential for r in rows)
        sub = filter(r -> r.potential == pot, rows)
        nus = sort!(unique(r.nu for r in sub))
        plt = plot(; xlabel = "perturbation order N",
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
    isempty(rows) && return nothing
    plt = plot(; xlabel = "perturbation order N",
               ylabel = "Mathematica / Julia (median time ratio)",
               yscale = :log10,
               title  = "Speedup factor — exact rational arithmetic",
               legend = :topright,
               size   = (800, 500))
    for pot in unique(r.potential for r in rows)
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
    n_total = length(rows)
    n_match = count(r -> r.matched, rows)

    io = IOBuffer()
    println(io, "# BenderWu — Julia vs Mathematica benchmarks\n")
    println(io, "Comparison of this Julia package against the reference Mathematica")
    println(io, "implementation `BenderWu.m` (arXiv:1608.08256). Both")
    println(io, "implementations compute the perturbative energy corrections ε_l for")
    println(io, "l = 0…N at fixed quantum number ν using the polynomial potentials")
    println(io, "shown below.\n")

    println(io, "Only **exact-rational arithmetic** is benchmarked. A Float64 vs")
    println(io, "MachinePrecision comparison would not be apples-to-apples:")
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

    println(io, "Julia timings come from `BenchmarkTools.@benchmark` (median over")
    println(io, "many samples, with a fresh `Potential` per sample so caches are")
    println(io, "cold). Mathematica timings come from `RepeatedTiming` which")
    println(io, "averages an automatically chosen number of repetitions.\n")

    println(io, "## Validation\n")
    println(io, "Both implementations were checked to agree on ε_l (even orders)")
    println(io, "before timings were recorded — bit-for-bit equality on rationals.\n")
    println(io, "**$(n_match)/$(n_total)** cases match exactly.\n")

    if n_match < n_total
        println(io, "Mismatched cells:\n")
        println(io, "| potential | ν | N |")
        println(io, "|---|---|---|")
        for r in rows
            r.matched && continue
            println(io, "| $(r.potential) | $(r.nu) | $(r.N) |")
        end
        println(io)
    end

    println(io, "## Results\n")
    println(io, "| potential | ν | N | Julia | Mathematica | speedup |")
    println(io, "|---|---|---|---:|---:|---:|")
    for r in sort(rows, by = r -> (r.potential, r.nu, r.N))
        println(io, "| $(r.potential) | $(r.nu) | $(r.N) | ",
                fmt_ms(r.jl_ms), " ms | ",
                fmt_ms(r.mma_ms), " ms | ",
                fmt_ratio(r.ratio), " |")
    end
    println(io)

    println(io, "### Per-potential timings\n")
    for pot in sort!(unique(r.potential for r in rows))
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
    println("Joined ", length(rows), " cases shared between Julia and Mathematica.")
    n_match = count(r -> r.matched, rows)
    println("  $n_match/$(length(rows)) match exactly.")
    for r in rows
        r.matched && continue
        @printf("  MISMATCH: %s ν=%d N=%d\n", r.potential, r.nu, r.N)
    end

    plot_per_potential(rows)
    plot_speedup(rows)

    jl_machine  = Dict(string(k) => v for (k, v) in pairs(jl.machine))
    mma_machine = Dict(string(k) => v for (k, v) in pairs(mma.machine))

    outpath = joinpath(ROOT, "BENCHMARKS.md")
    write_report(rows, jl_machine, mma_machine; outpath)
    println("Wrote ", outpath)
end

main()
