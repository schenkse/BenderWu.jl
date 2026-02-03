# Helper functions
# [1608.08256]

# Maximum number K_l
@memoize function max_k(ν::Int, l::Int, vcoeffs)
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

@memoize function A_kl(ν::Int, k::Int, l::Int, vcoeffs)
    ω = sqrt(2 * vcoeffs[1])

    # Minimum value for k or l
    if k < 0 || l < 0 return zero(ω) end

    # Maximum value for k
    if k > ν && iszero(l) return zero(ω) end
    if k == ν && iszero(l) return one(ω) end
    if k == ν && l > 0 return zero(ω) end
    if k > max_k(ν, l, vcoeffs) return zero(ω) end
    
    Akl = (k+2) * (k+1) * A_kl(ν, k+2, l, vcoeffs)
    if k > ν && l > 0
        # Terminate sum for a finite number of terms in the potential
        for n=1:l
            if n != l
                Akl += 2 * ε_l(ν, n, vcoeffs) * A_kl(ν, k, l-n, vcoeffs)
            end
            if n+1 > length(vcoeffs) continue end
            if iszero(vcoeffs[n+1]) continue end
            Akl += -2 * vcoeffs[n+1] * A_kl(ν, k-n-2, l-n, vcoeffs)
        end
    else
        # Terminate sum for a finite number of terms in the potential
        for n=1:l
            Akl += 2 * ε_l(ν, n, vcoeffs) * A_kl(ν, k, l-n, vcoeffs)
            if n+1 > length(vcoeffs) continue end
            if iszero(vcoeffs[n+1]) continue end
            Akl += -2 * vcoeffs[n+1] * A_kl(ν, k-n-2, l-n, vcoeffs)
        end
    end
    return Akl / (2 * ω * (k - ν))
end

@memoize function ε_l(ν::Int, l::Int, vcoeffs)
    ω = sqrt(2 * vcoeffs[1])
    if isodd(l) return zero(ω) end
    if iszero(l) return ω * (ν + 1/2) end
    ε = -(ν+2) * (ν+1) / 2 * A_kl(ν, ν+2, l, vcoeffs)
    # Terminate sum for a finite number of terms in the potential
    for n=1:l
        if n+1 > length(vcoeffs) continue end
        if iszero(vcoeffs[n+1]) continue end
        ε += vcoeffs[n+1] * A_kl(ν, ν-n-2, l-n, vcoeffs)
    end
    return ε
end

# Here goes some iterative solution to the same problem
function initialize_Akl_eps(ν::Int, l::Int, vcoeffs)
    kmax = max_k(ν, l, vcoeffs)
    otype = typeof(vcoeffs[1])
    return zeros(otype, kmax+3, l+1), zeros(otype, l+1)
end

function fill_Akl!(Akl, ε, ν::Int, maxorder::Int, vcoeffs)
    ω = sqrt(2 * vcoeffs[1])
    # Be careful with indexing here
    Akl[ν+1, 1] = one(ω)
    ε[1] = ω * (ν + 1/2)
    for l=0:maxorder
        kmax = max_k(ν, l, vcoeffs)
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
function find_epoly(order::Int, vcoeffs)
    if isodd(order)
        return zeros(typeof(vcoeffs[1]), floor(Int, order/2) + 2)
    end
    # At order l we need to compute l+2 terms in total
    l = Int(order/2)
    ε_n = [ε_l(n, order, vcoeffs) for n=0:l+1]
    N_mat = [oftype(vcoeffs[1], n^j) for n=0:l+1, j=0:l+1]
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