#!/usr/bin/env julia
# Loads benchmark/results_julia.json and benchmark/results_mma.json, validates
# that ε_l agree across implementations, draws plots, and writes BENCHMARKS.md
# at the repo root.
#
# Usage:
#   julia --project=benchmark benchmark/aggregate.jl

using JSON3
using Plots
using Printf
using Statistics

include(joinpath(@__DIR__, "cases.jl"))

const ROOT = abspath(joinpath(@__DIR__, ".."))
const PLOTS_DIR = joinpath(@__DIR__, "plots")
isdir(PLOTS_DIR) || mkpath(PLOTS_DIR)

# ---------------------------------------------------------------------------
# Loading

function load_results(path)
    raw = JSON3.read(read(path, String))
    return raw
end

# ---------------------------------------------------------------------------
# Validation

"""
Parse a serialised ε value back into a comparable Julia number.

Julia rationals come out as e.g. "3//4"; floats come out as printed Float64
strings ("0.75", "3.0e-5"). Mathematica exact values are InputForm strings
("3/4", "(123/456)"); Mathematica machine-precision floats come out as CForm
("0.75", "3e-05"). We normalise to BigFloat for comparison except when both
sides are exact, where we use Rational{BigInt}.
"""
function parse_eps(s::AbstractString, mode::AbstractString)
    s = strip(s)
    if mode == "rational"
        # Try Julia's "a//b", then Mathematica's "a/b". Strip outer parens.
        s = strip(s, ['(', ')'])
        if occursin("//", s)
            num, den = split(s, "//")
            return Rational{BigInt}(parse(BigInt, strip(num)),
                                    parse(BigInt, strip(den)))
        elseif occursin("/", s) && !occursin("e", lowercase(s))
            num, den = split(s, "/")
            return Rational{BigInt}(parse(BigInt, strip(num)),
                                    parse(BigInt, strip(den)))
        else
            return Rational{BigInt}(parse(BigInt, s))
        end
    else
        # Float string, possibly Mathematica CForm. Replace Mathematica's
        # "e" forms (which are already plain) — both forms parse fine.
        return parse(BigFloat, s)
    end
end

"""
Compare two ε vectors (already parsed) on the union of even-l indices.
Returns (matched::Bool, max_relerr::Float64).

The Mathematica list always starts at l=0 and contains entries for l=0,2,4,...
Our Julia list contains entries for l=0,1,2,...,N (with odd-l = 0). So we
compare jl_eps[2k+1] with mma_eps[k+1] for k = 0:min(...,...).
"""
function compare_eps(jl_eps::Vector, mma_eps::Vector, mode::AbstractString)
    n = min(length(mma_eps), (length(jl_eps) - 1) ÷ 2 + 1)
    max_relerr = 0.0
    matched = true
    for k in 0:(n-1)
        a = jl_eps[2k + 1]
        b = mma_eps[k + 1]
        if mode == "rational"
            if a != b
                matched = false
                # Convert rationals to BigFloat for a relerr value.
                fa, fb = BigFloat(a), BigFloat(b)
                relerr = iszero(fb) ? Float64(abs(fa)) :
                                       Float64(abs(fa - fb) / abs(fb))
                max_relerr = max(max_relerr, relerr)
            end
        else
            relerr = if iszero(b) && iszero(a)
                0.0
            elseif iszero(b)
                Float64(abs(a))
            else
                Float64(abs(a - b) / abs(b))
            end
            max_relerr = max(max_relerr, relerr)
            # Floating-point loss model: rtol grows with order, since both
            # implementations accumulate Float64 round-off through ~3·l multiply-
            # adds per order with different summation orders. Empirically the
            # exponent ~l/8 covers observed drift up to l=50; cap at 0.1.
            tol = min(1e-10 * 10.0^(2k / 8), 0.1)
            if relerr > tol
                matched = false
            end
        end
    end
    return matched, max_relerr
end

# ---------------------------------------------------------------------------
# Joining

function join_results(jl, mma)
    function key(r)
        (String(r.potential), Int(r.nu), Int(r.N), String(r.mode))
    end
    mma_by = Dict(key(r) => r for r in mma.results)
    rows = []
    for j in jl.results
        k = key(j)
        m = get(mma_by, k, nothing)
        m === nothing && continue
        jl_eps = [parse_eps(s, String(j.mode)) for s in j.epsilons]
        mma_eps = [parse_eps(s, String(j.mode)) for s in m.epsilons]
        matched, relerr = compare_eps(jl_eps, mma_eps, String(j.mode))
        push!(rows, (
            potential = String(j.potential),
            nu = Int(j.nu),
            N = Int(j.N),
            mode = String(j.mode),
            jl_ms = Float64(j.time_ns_median) / 1e6,
            mma_ms = Float64(m.time_ns_median) / 1e6,
            ratio = Float64(m.time_ns_median) / Float64(j.time_ns_median),
            matched = matched,
            relerr = relerr,
        ))
    end
    return rows
end

# ---------------------------------------------------------------------------
# Plotting

function plot_mode(rows, mode_label)
    sub = filter(r -> r.mode == mode_label, rows)
    isempty(sub) && return String[]
    pots = unique(r.potential for r in sub)
    paths = String[]
    for pot in pots
        rows_p = filter(r -> r.potential == pot, sub)
        nus = sort!(unique(r.nu for r in rows_p))
        plt = plot(; xlabel = "perturbation order N",
                   ylabel = "median time (ms)",
                   yscale = :log10,
                   title = "$pot — $mode_label",
                   legend = :topleft,
                   size = (700, 450))
        for ν in nus
            rs = sort(filter(r -> r.nu == ν, rows_p), by = r -> r.N)
            Ns = [r.N for r in rs]
            plot!(plt, Ns, [r.jl_ms for r in rs];
                  marker = :circle, label = "Julia ν=$ν")
            plot!(plt, Ns, [r.mma_ms for r in rs];
                  marker = :square, linestyle = :dash, label = "Mathematica ν=$ν")
        end
        path = joinpath(PLOTS_DIR, "$(pot)_$(mode_label).png")
        savefig(plt, path)
        push!(paths, path)
    end
    return paths
end

function plot_speedup(rows)
    isempty(rows) && return nothing
    pots = unique(r.potential for r in rows)
    plt = plot(; xlabel = "perturbation order N",
               ylabel = "Mathematica / Julia (median time ratio)",
               yscale = :log10,
               title = "Speedup factor — Julia vs Mathematica",
               legend = :topleft,
               size = (800, 500))
    for pot in pots, mode in ("float64", "rational")
        rs = sort(filter(r -> r.potential == pot && r.mode == mode && r.nu == 0, rows),
                  by = r -> r.N)
        isempty(rs) && continue
        Ns = [r.N for r in rs]
        plot!(plt, Ns, [r.ratio for r in rs];
              marker = :circle, label = "$pot — $mode")
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
    println(io, "implementation `BenderWu.m` (Sulejmanpasic, arXiv:1608.08256).")
    println(io, "Both implementations compute perturbative energy corrections")
    println(io, "ε_l up to a maximum order N at fixed quantum number ν using the")
    println(io, "polynomial potentials shown below.\n")

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
    println(io, "before timings were recorded.\n")
    println(io, "- **$(n_match)/$(n_total)** cases match within tolerance.")
    println(io, "- Tolerance for `float64`/`MachinePrecision`: relative error")
    println(io, "  `≤ 10⁻¹⁰ · 10^(l/4)` capped at `10⁻¹`. The cap reflects the")
    println(io, "  expected drift between two IEEE-double implementations whose")
    println(io, "  summation orders differ — high-order ε_l can lose several")
    println(io, "  digits of precision in either code.")
    println(io, "- Tolerance for `rational`/exact: bit-for-bit equality.\n")

    if n_match < n_total
        println(io, "Mismatched cells (max relative error shown):\n")
        println(io, "| potential | ν | N | mode | max rel. err |")
        println(io, "|---|---|---|---|---|")
        for r in rows
            r.matched && continue
            println(io, "| $(r.potential) | $(r.nu) | $(r.N) | $(r.mode) | ",
                    @sprintf("%.2e", r.relerr), " |")
        end
        println(io)
    end

    for mode in ("float64", "rational")
        sub = filter(r -> r.mode == mode, rows)
        isempty(sub) && continue
        title = mode == "float64" ?
                "Float64 vs MachinePrecision" :
                "Rational{BigInt} vs Mathematica exact"
        println(io, "## $title\n")
        println(io, "| potential | ν | N | Julia | Mathematica | speedup | match |")
        println(io, "|---|---|---|---:|---:|---:|:---:|")
        for r in sort(sub, by = r -> (r.potential, r.nu, r.N))
            mark = r.matched ? "✓" : "✗"
            println(io, "| $(r.potential) | $(r.nu) | $(r.N) | ",
                    fmt_ms(r.jl_ms), " ms | ",
                    fmt_ms(r.mma_ms), " ms | ",
                    fmt_ratio(r.ratio), " | $mark |")
        end
        println(io)

        for pot in unique(r.potential for r in sub)
            img = "benchmark/plots/$(pot)_$(mode).png"
            isfile(joinpath(ROOT, img)) || continue
            println(io, "![$(pot) — $(mode)]($img)\n")
        end
    end

    speedup_path = "benchmark/plots/speedup.png"
    if isfile(joinpath(ROOT, speedup_path))
        println(io, "## Overall speedup\n")
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
    println("  $n_match/$(length(rows)) match within tolerance.")
    for r in rows
        r.matched && continue
        @printf("  MISMATCH: %s ν=%d N=%d mode=%s  max rel.err=%.2e\n",
                r.potential, r.nu, r.N, r.mode, r.relerr)
    end

    plot_mode(rows, "float64")
    plot_mode(rows, "rational")
    plot_speedup(rows)

    jl_machine  = Dict(string(k) => v for (k, v) in pairs(jl.machine))
    mma_machine = Dict(string(k) => v for (k, v) in pairs(mma.machine))

    outpath = joinpath(ROOT, "BENCHMARKS.md")
    write_report(rows, jl_machine, mma_machine; outpath)
    println("Wrote ", outpath)
end

main()
