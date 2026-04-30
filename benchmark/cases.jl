# Shared definition of the benchmark matrix used by run_julia.jl and aggregate.jl.
# The Mathematica driver mirrors this list in run_mathematica.wls.

# A potential is described by its Julia coefficient vector (vcoeffs[n] is the
# coefficient of x^(n+1)) and a Mathematica string representation of V(x).
const POTENTIALS = [
    (name = "quartic",      vcoeffs = [1//2, 0//1, 1//1],                              mma = "x^2/2 + x^4"),
    (name = "sextic",       vcoeffs = [1//2, 0//1, 0//1, 0//1, 1//1],                  mma = "x^2/2 + x^6"),
    (name = "octic",        vcoeffs = [1//2, 0//1, 0//1, 0//1, 0//1, 0//1, 1//1],      mma = "x^2/2 + x^8"),
    (name = "mixed_parity", vcoeffs = [1//2, 1//1, 1//1],                              mma = "x^2/2 + x^3 + x^4"),
]

const NUS = [0, 1, 5]

# Order sweeps. Float64 holds up well; rational/exact balloons the integers.
const ORDERS_FLOAT   = [10, 20, 30, 40, 50]
const ORDERS_RATIONAL = [5, 10, 15, 20, 25]

const QUICK_NUS = [0]
const QUICK_ORDERS_FLOAT = [10]
const QUICK_ORDERS_RATIONAL = [5]

"""
    cases(; quick=false)

Return a vector of named tuples describing every (potential, ν, N, mode) cell
in the benchmark matrix. `mode` is `:float64` or `:rational`.
"""
function cases(; quick::Bool = false)
    nus = quick ? QUICK_NUS : NUS
    orders_f = quick ? QUICK_ORDERS_FLOAT : ORDERS_FLOAT
    orders_r = quick ? QUICK_ORDERS_RATIONAL : ORDERS_RATIONAL
    out = Vector{NamedTuple}()
    for p in POTENTIALS, ν in nus
        for N in orders_f
            push!(out, (potential = p.name, vcoeffs = p.vcoeffs, mma = p.mma,
                        nu = ν, N = N, mode = :float64))
        end
        for N in orders_r
            push!(out, (potential = p.name, vcoeffs = p.vcoeffs, mma = p.mma,
                        nu = ν, N = N, mode = :rational))
        end
    end
    return out
end

case_key(c) = (c.potential, c.nu, c.N, String(c.mode))
