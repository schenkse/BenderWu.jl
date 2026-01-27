# Helper functions
# [1608.08256]

# Maximum number K_l
function max_kl(ν::Int, l::Int, vcoeffs)
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

function A_kl(ν::Int, k::Int, l::Int, vcoeffs)
    ω = sqrt(2 * vcoeffs[1])
    # Maximum value for k
    if k > max_kl(ν, l, vcoeffs) return 0.0 end
    if k > ν && iszero(l) return 0.0 end
    if k == ν && iszero(l) return 1.0 end
    if k == ν && l > 0 return 0.0 end

    # Minimum value for k or l
    if k < 0 || l < 0 return 0.0 end
    
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
    end
    if k < ν && l > 0
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

function ε_l(ν::Int, l::Int, vcoeffs)
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