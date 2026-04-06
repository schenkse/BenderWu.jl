module BenderWu

# Bender-Wu method for perturbative energy levels of polynomial potentials
# Reference: arXiv:1608.08256

export Potential
export max_k, A_kl, ε_l
export initialize_Akl_eps, fill_Akl!
export find_epoly, find_epoly_derivative, evaluate_epoly

"""
    Potential(vcoeffs)

Represents a polynomial potential with coefficients `vcoeffs`, where `vcoeffs[n]`
is the coefficient of x^(n-1). Carries its own memoization caches, which are
GC-managed — create one instance per potential and reuse it across calls.

# Example
```julia
pot = Potential([0.5, 0.0, 1.0])  # V(x) = x²/2 + x⁴
epoly = find_epoly(2, pot)
```
"""
struct Potential{T}
    vcoeffs::Vector{T}
    _max_k_cache::Dict{Tuple{Int,Int}, Int}
    _Akl_cache::Dict{Tuple{Int,Int,Int}, T}
    _εl_cache::Dict{Tuple{Int,Int}, T}
end

Potential(vcoeffs::AbstractVector{T}) where T = Potential(
    collect(T, vcoeffs),
    Dict{Tuple{Int,Int}, Int}(),
    Dict{Tuple{Int,Int,Int}, T}(),
    Dict{Tuple{Int,Int}, T}()
)

# Maximum number K_l
function max_k(pot::Potential, ν::Int, l::Int)
    get!(pot._max_k_cache, (ν, l)) do
        vcoeffs = pot.vcoeffs
        L = findfirst(!iszero, vcoeffs[2:end])
        if 0 <= l < L
            if iszero(l)
                return ν
            else
                return 0
            end
        elseif l >= L
            return ν + (L + 2) * floor(Int, l / L) + l % L
        end
    end
end

function A_kl(pot::Potential, ν::Int, k::Int, l::Int)
    T = eltype(pot.vcoeffs)

    # Cheap boundary cases — not worth caching
    if k < 0 || l < 0 return zero(T) end
    if k > ν && iszero(l) return zero(T) end
    if k == ν && iszero(l) return one(T) end
    if k == ν && l > 0 return zero(T) end
    if k > max_k(pot, ν, l) return zero(T) end

    get!(pot._Akl_cache, (ν, k, l)) do
        vcoeffs = pot.vcoeffs
        ω = sqrt(2 * vcoeffs[1])
        Akl = (k+2) * (k+1) * A_kl(pot, ν, k+2, l)
        if k > ν && l > 0
            # Terminate sum for a finite number of terms in the potential
            for n=1:l
                if n != l
                    Akl += 2 * ε_l(pot, ν, n) * A_kl(pot, ν, k, l-n)
                end
                if n+1 > length(vcoeffs) continue end
                if iszero(vcoeffs[n+1]) continue end
                Akl += -2 * vcoeffs[n+1] * A_kl(pot, ν, k-n-2, l-n)
            end
        else
            # Terminate sum for a finite number of terms in the potential
            for n=1:l
                Akl += 2 * ε_l(pot, ν, n) * A_kl(pot, ν, k, l-n)
                if n+1 > length(vcoeffs) continue end
                if iszero(vcoeffs[n+1]) continue end
                Akl += -2 * vcoeffs[n+1] * A_kl(pot, ν, k-n-2, l-n)
            end
        end
        return Akl / (2 * ω * (k - ν))
    end
end

function ε_l(pot::Potential, ν::Int, l::Int)
    T = eltype(pot.vcoeffs)
    ω = sqrt(2 * pot.vcoeffs[1])

    # Cheap boundary cases — not worth caching
    if isodd(l) return zero(T) end
    if iszero(l) return ω * (ν + one(T)/2) end

    get!(pot._εl_cache, (ν, l)) do
        vcoeffs = pot.vcoeffs
        ε = -(ν+2) * (ν+1) / 2 * A_kl(pot, ν, ν+2, l)
        # Terminate sum for a finite number of terms in the potential
        for n=1:l
            if n+1 > length(vcoeffs) continue end
            if iszero(vcoeffs[n+1]) continue end
            ε += vcoeffs[n+1] * A_kl(pot, ν, ν-n-2, l-n)
        end
        return ε
    end
end

# Iterative solution
function initialize_Akl_eps(pot::Potential, ν::Int, l::Int)
    kmax = max_k(pot, ν, l)
    T = eltype(pot.vcoeffs)
    return zeros(T, kmax+3, l+1), zeros(T, l+1)
end

function fill_Akl!(Akl, ε, pot::Potential, ν::Int, maxorder::Int)
    vcoeffs = pot.vcoeffs
    ω = sqrt(2 * vcoeffs[1])
    # Be careful with indexing here
    Akl[ν+1, 1] = one(ω)
    ε[1] = ω * (ν + 1/2)
    for l=0:maxorder
        kmax = max_k(pot, ν, l)
        # Step 1
        if l > 0
            for k=kmax:-1:ν+1
                Akl[k+1, l+1] += (k+2) * (k+1) * Akl[k+2+1, l+1]
                for n=1:l-1
                    Akl[k+1, l+1] += 2 * ε[n+1] * Akl[k+1, l-n+1]
                end
                for n=1:l
                    # Check index bounds
                    if n+1 > length(vcoeffs) || k-n-2 < 0 continue end
                    Akl[k+1, l+1] += -2 * vcoeffs[n+1] * Akl[k-n-2+1, l-n+1]
                end
                Akl[k+1, l+1] /= 2 * ω * (k - ν)
            end
        end
        # Step 2
        if l > 0
            ε[l+1] += -(ν+2) * (ν+1) / 2 * Akl[ν+2+1, l+1]
            for n=1:l
                # Check index bounds
                if n+1 > length(vcoeffs) || ν-n-2 < 0 continue end
                ε[l+1] += vcoeffs[n+1] * Akl[ν-n-2+1, l-n+1]
            end
        end
        # Step 3
        for k=ν-1:-1:0
            Akl[k+1, l+1] += (k+2) * (k+1) * Akl[k+2+1, l+1]
            for n=1:l
                Akl[k+1, l+1] += 2 * ε[n+1] * Akl[k+1, l-n+1]
                # Check index bounds
                if n+1 > length(vcoeffs) || k-n-2 < 0 continue end
                Akl[k+1, l+1] += -2 * vcoeffs[n+1] * Akl[k-n-2+1, l-n+1]
            end
            Akl[k+1, l+1] /= 2 * ω * (k - ν)
        end
    end
    nothing
end

# Fit polynomial to energy coefficients
function find_epoly(order::Int, pot::Potential)
    if isodd(order)
        return zeros(eltype(pot.vcoeffs), floor(Int, order/2) + 2)
    end
    # At order l we need to compute l+2 terms in total
    l = Int(order/2)
    ε_n = [ε_l(pot, n, order) for n=0:l+1]
    N_mat = [oftype(pot.vcoeffs[1], n^j) for n=0:l+1, j=0:l+1]
    return N_mat \ ε_n
end

function find_epoly_derivative(epoly)
    otype = typeof(epoly[1])
    ds = Array{otype}(undef, length(epoly)-1)
    for (n, ε) in enumerate(epoly)
        if isone(n) continue end
        ds[n-1] = n >= 20 ? factorial(big(n-1)) * ε : factorial(n-1) * ε
    end
    return ds
end

function evaluate_epoly(n::Int, epoly)
    res = zero(epoly[1])
    for (k, ε) in enumerate(epoly)
        res += ε * n^(k-1)
    end
    return res
end

end # module BenderWu
