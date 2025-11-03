import Base: *

struct Symmetry
    sx::Int64   # should have inner constructor to ensure s,sx,sy,sz = ±1
    sy::Int64
    sz::Int64
    ax::Rational{Int}
    az::Rational{Int}
    s::Int64
end

Symmetry() = Symmetry(1, 1, 1, 0//1, 0//1, 1) # the identity
Symmetry(sx,sy,sz) = Symmetry(sx, sy, sz, 0//1, 0//1, 1)
Symmetry(sx,sy,sz,ax,az) = Symmetry(sx, sy, sz, ax, az, 1)

*(σ1::Symmetry, σ2::Symmetry) = Symmetry(σ1.sx*σ2.sx, σ1.sy*σ2.sy, σ1.sz*σ2.sz, 
                                         (σ1.ax + σ2.ax) % 1, (σ1.az + σ2.az) % 1, σ1.s*σ2.s)
    

"""
   symmetric(ijkl::Vector{Int}, σ::Symmetry)

Return true if σ Ψᵢⱼₖₗ = Ψᵢⱼₖₗ 
"""
function symmetric(ijkl::Vector{Int}, σ::Symmetry)
    
    rtn = σ.s
    rtn *= σ.sx == -1 ? xreflection(ijkl) : 1
    rtn *= σ.sy == -1 ? yreflection(ijkl) : 1
    rtn *= σ.sz == -1 ? zreflection(ijkl) : 1

    if σ.ax == 1//2
        rtn *= xtranslationLx2(ijkl)
    elseif σ.ax != 0//2
        error("symmetries only implemented for half-box shifts")
    end

    if σ.az == 1//2
        rtn *= ztranslationLz2(ijkl)
    elseif σ.az != 0//2
        error("symmetries only implemented for half-box shifts")
    end

    return rtn == 1
end


"""
   symmetric(ijkl::Vector{Int}, σ::Vector{Symmetry})

Return true if σ[n] Ψᵢⱼₖₗ = Ψᵢⱼₖₗ n=1:end
"""
function symmetric(ijkl::Vector{Int}, σ::Vector{Symmetry})
    for n=1:length(σ)
        if !symmetric(ijkl, σ[n])
            return false
        end
    end
    true
end


