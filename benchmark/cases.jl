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

# A potential is described by its Julia coefficient vector (vcoeffs[n] is the
# coefficient of x^(n+1)) and a Mathematica string representation of V(x).
const POTENTIALS = [
    (name = "quartic",      vcoeffs = [1//2, 0//1, 1//1],                              mma = "x^2/2 + x^4"),
    (name = "sextic",       vcoeffs = [1//2, 0//1, 0//1, 0//1, 1//1],                  mma = "x^2/2 + x^6"),
    (name = "octic",        vcoeffs = [1//2, 0//1, 0//1, 0//1, 0//1, 0//1, 1//1],      mma = "x^2/2 + x^8"),
    (name = "mixed_parity", vcoeffs = [1//2, 1//1, 1//1],                              mma = "x^2/2 + x^3 + x^4"),
]

const NUS = [0, 1, 5]

# Rational/exact arithmetic: integer growth is super-exponential in N, so the
# range stays modest.
const ORDERS = [5, 10, 15, 20, 25]

const QUICK_NUS = [0]
const QUICK_ORDERS = [5]

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
