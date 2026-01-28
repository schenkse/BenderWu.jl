# Helper functions
# [1608.08256]

# Maximum number K_l
@memoize Dict function max_k(ν::Int, l::Int, vcoeffs)
    L = findfirst(!iszero, vcoeffs[2:end]) - 1
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

@memoize Dict function A_kl(ν::Int, k::Int, l::Int, vcoeffs)
    ω = sqrt(2 * vcoeffs[1])

    # Minimum value for k or l
    if k < 0 || l < 0 return 0.0 end

    # Maximum value for k
    if k > max_k(ν, l, vcoeffs) return 0.0 end
    if k > ν && iszero(l) return 0.0 end
    if k == ν && iszero(l) return 1.0 end
    if k == ν && l > 0 return 0.0 end
    
    if k > ν && l > 0
        Akl = (k+2) * (k+1) * A_kl(ν, k+2, l, vcoeffs)
        # Terminate sum for a finite number of terms in the potential
        lmin = min(l, length(vcoeffs)-1)
        for n=1:lmin
            if iszero(vcoeffs[n+1]) continue end
            Akl += -2 * vcoeffs[n+1] * A_kl(ν, k-n-2, l-n, vcoeffs)
        end
        for n=1:l-1
            Akl += 2 * ε_l(ν, n, vcoeffs) * A_kl(ν, k, l-n, vcoeffs)
        end
        
        return Akl / (2 * ω * (k - ν))
    else
        Akl = (k+2) * (k+1) * A_kl(ν, k+2, l, vcoeffs)
        # Terminate sum for a finite number of terms in the potential
        lmin = min(l, length(vcoeffs)-1)
        for n=1:lmin
            if !iszero(vcoeffs[n+1])
                Akl += -2 * vcoeffs[n+1] * A_kl(ν, k-n-2, l-n, vcoeffs)
            end
            Akl += 2 * ε_l(ν, n, vcoeffs) * A_kl(ν, k, l-n, vcoeffs)
        end
        return Akl / (2 * ω * (k - ν))
    end
end

@memoize Dict function ε_l(ν::Int, l::Int, vcoeffs)
    ω = sqrt(2 * vcoeffs[1])
    if iszero(l) return ω * (ν + 1/2) end
    ε = -(ν+2) * (ν+1) / 2 * A_kl(ν, ν+2, l, vcoeffs)
    # Terminate sum for a finite number of terms in the potential
    lmin = min(l, length(vcoeffs)-1)
    for n=1:lmin
        if iszero(vcoeffs[n+1]) continue end
        ε += vcoeffs[n+1] * A_kl(ν, ν-n-2, l-n, vcoeffs)
    end
    return ε
end

# Here goes some iterative solution to the same problem
function initialize_Akl(ν::Int, l::Int, vcoeffs)
    kmax = max_k(ν, l, vcoeffs)
    return zeros(kmax+3, l+1)
end
initialize_eps(ν::Int, l::Int) = zeros(l+1)

function fill_Akl!(Akl, ε, ν::Int, maxorder::Int, vcoeffs)
    ω = sqrt(2.0 * vcoeffs[1])
    # Be careful with indexing here
    Akl[ν+1, 1] = 1.0
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
                    #if iszero(vcoeffs[n+1]) continue end
                    Akl[k+1, l+1] += -2 * vcoeffs[n+1] * Akl[k-n-2+1, l-n+1]
                end
                Akl[k+1, l+1] /= 2 * ω * (k - ν)
            end
        end
        # Step 2
        if l > 0
            ε[l+1] += -(ν+2) * (ν+1) / 2 * Akl[ν+2+1, l+1]
            for n=1:l
                #if iszero(vcoeffs[n+1]) continue end
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
                #if iszero(vcoeffs[n+1]) continue end
                # Check index bounds
                if n+1 > length(vcoeffs) || k-n-2 < 0 continue end
                Akl[k+1, l+1] += -2 * vcoeffs[n+1] * Akl[k-n-2+1, l-n+1]
            end
            Akl[k+1, l+1] /= 2 * ω * (k - ν)
        end
    end
    nothing
end