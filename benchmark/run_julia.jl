#!/usr/bin/env julia
# Times the Julia BenderWu implementation across the shared benchmark matrix
# and writes results to benchmark/results_julia.json.
#
# Usage:
#   julia --project=benchmark benchmark/run_julia.jl [--quick]

using BenchmarkTools
using BenderWu
using JSON3

include(joinpath(@__DIR__, "cases.jl"))

const QUICK = "--quick" in ARGS

# BenchmarkTools defaults: keep total time per case bounded so the rational
# sweeps stay tractable. Each case still gets multiple samples for a stable
# median.
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5.0
BenchmarkTools.DEFAULT_PARAMETERS.samples = 50

"""
Build a Potential of the right element type and time the iterative
fill_Akl! call. Each sample uses a fresh Potential so caches are cold —
this matches what Mathematica does (every call recomputes from scratch).
"""
function bench_case(c)
    coeffs = c.mode === :float64 ? Float64.(c.vcoeffs) :
                                    Rational{BigInt}.(c.vcoeffs)
    ν = c.nu
    N = c.N

    # Reference computation: get the ε vector once, with full precision
    # preserved (rationals as exact). This is what we serialise for the
    # validation pass.
    pot_ref = Potential(coeffs)
    Akl_ref, ε_ref = initialize_Akl_eps(pot_ref, ν, N)
    fill_Akl!(Akl_ref, ε_ref, pot_ref, ν, N)

    # Timed call: rebuild Potential and arrays every sample.
    b = @benchmark begin
        pot = Potential($coeffs)
        Akl, ε = initialize_Akl_eps(pot, $ν, $N)
        fill_Akl!(Akl, ε, pot, $ν, $N)
    end

    return (
        potential = c.potential,
        nu = ν,
        N = N,
        mode = String(c.mode),
        time_ns_median = median(b.times),
        time_ns_min = minimum(b.times),
        samples = length(b.times),
        allocs = b.allocs,
        memory_bytes = b.memory,
        epsilons = [string(e) for e in ε_ref],
    )
end

function main()
    cs = cases(; quick = QUICK)
    println("Running ", length(cs), " Julia cases", QUICK ? " (quick mode)" : "", "...")
    results = []
    for (i, c) in enumerate(cs)
        print("[", i, "/", length(cs), "] ", c.potential, " ν=", c.nu,
              " N=", c.N, " mode=", c.mode, " ... ")
        flush(stdout)
        r = bench_case(c)
        push!(results, r)
        println(round(r.time_ns_median / 1e6, digits=3), " ms (median)")
    end

    machine = Dict(
        "julia_version" => string(VERSION),
        "cpu" => Sys.cpu_info()[1].model,
        "ncores" => Sys.CPU_THREADS,
        "os" => string(Sys.KERNEL),
        "arch" => string(Sys.ARCH),
    )

    out = Dict("machine" => machine, "results" => results)
    open(joinpath(@__DIR__, "results_julia.json"), "w") do io
        JSON3.pretty(io, out)
    end
    println("Wrote ", joinpath("benchmark", "results_julia.json"))
end

main()
