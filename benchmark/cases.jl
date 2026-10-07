# Shared definition of the benchmark matrix used by run_julia.jl and aggregate.jl.
# The Mathematica driver mirrors this list in run_mathematica.wls.
#
# We compare only exact-rational arithmetic (Julia Rational{BigInt} vs
# Mathematica's exact symbolic arithmetic). A Float64 vs MachinePrecision
# comparison would not be apples-to-apples: Mathematica's BenderWu evaluates
# the recursion through its symbolic term-rewriting pipeline regardless of
# the coefficient precision, while Julia compiles to tight native loops.
# That comparison would mostly measure symbolic-evaluator overhead, not
# algorithmic efficiency, so we omit it.

#
# N is the highest power of g in the expansion E = Σ ε_l g^l. Julia computes
# l = 0:N, while Mathematica's BenderWu[V, x, ν, n] counts powers of g², so it
# is called with n = N/2. Every N must therefore be even.

# A potential is described by its Julia coefficient vector (vcoeffs[n] is the
# coefficient of x^(n+1)), a Mathematica string representation of V(x), and a
# display form for the report. The order here is the order in BENCHMARKS.md.
const POTENTIALS = [
    (name = "quartic",      vcoeffs = [1//2, 0//1, 1//1],                         mma = "x^2/2 + x^4",       display = "x²/2 + x⁴"),
    (name = "mixed_parity", vcoeffs = [1//2, 1//1, 1//1],                         mma = "x^2/2 + x^3 + x^4", display = "x²/2 + x³ + x⁴"),
    (name = "sextic",       vcoeffs = [1//2, 0//1, 0//1, 0//1, 1//1],             mma = "x^2/2 + x^6",       display = "x²/2 + x⁶"),
    (name = "octic",        vcoeffs = [1//2, 0//1, 0//1, 0//1, 0//1, 0//1, 1//1], mma = "x^2/2 + x^8",       display = "x²/2 + x⁸"),
]

const NUS = [0, 1, 5]

# Rational/exact arithmetic: integer growth is super-exponential in N, so the
# range stays modest.
const ORDERS = [10, 20, 30, 40, 50]

const QUICK_NUS = [0]
const QUICK_ORDERS = [10]

@assert all(iseven, ORDERS) && all(iseven, QUICK_ORDERS)

"""
    cases(; quick=false)

Return a vector of named tuples describing every (potential, ν, N) cell in the
benchmark matrix. All cases run in exact-rational mode.
"""
function cases(; quick::Bool = false)
    nus    = quick ? QUICK_NUS    : NUS
    orders = quick ? QUICK_ORDERS : ORDERS
    out = Vector{NamedTuple}()
    for p in POTENTIALS, ν in nus, N in orders
        push!(out, (potential = p.name, vcoeffs = p.vcoeffs, mma = p.mma,
                    nu = ν, N = N))
    end
    return out
end
