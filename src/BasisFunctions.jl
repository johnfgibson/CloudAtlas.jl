# divide(x::Int, y::Int) = x//y
# divide(x::Int, y::Rational) = x//y
# divide(x::Int, y::AbstractFloatRational) = x//y

# divide(x::Rational, y::Rational) = x//y
# divide(x::Rational, y::Int) = x//y
# divide(x::AbstractFloat, y::AbstractFloat) = x/y
# divide(x::AbstractFloat, y::AbstractFloat) = x/y

const T_ELEMENT = Float64

parity(k) = k%2 == 0 ? 1 : -1    # 1 for k even, -1 for k odd
square(x) = x*x

function polyparity(p::Polynomial)
    parity = 0     # none
    if sum(abs.(p.coeffs[1:2:end])) == 0
        parity = -1 # even
    elseif sum(abs.(p.coeffs[2:2:end])) == 0
        parity = 1  # odd
    end
    parity
end

"""
    FourierMode{T<:Real}

FourierMode represents a cos, sin or const function of form c * {cos(ajx) | 1 | sin(ajx)}.
Which of cos(ajx), 1, sin(ajx) a given FourierMode represents is determined by the index j
 
            cos αjx,  j<0                      cos γlz,  l<0
 E_j(x) = c 1,        j=0    and    E_l(z) = c 1,        l=0
            sin αjx,  j>0                      sin γlz,  l>0

The x and z are implicit; FourierModes are functions, don't know names of their operands
"""
struct FourierMode{T<:Real}
    coeff::T          # multiplicative coefficient
    waveindex::Int    # j of cos(ajx), and j<0 => cos(ajx), j=0 => 1, j>0 => sin(ajx)
    wavenumber::T     # alpha of cos(ajx), sin(ajx) or gamma of cos(gjx), sin(gjx)
    function FourierMode(coeff::TC, waveindex::Integer, wavenumber::TW) where {TC<:Real, TW<:Real}
        T = promote_type(TC, TW)
        return new{T}(
            T(coeff),
            waveindex,
            T(wavenumber)
        )
    end
end

# needed operations
#  zero?
#  isorthogonal
#  compatible
#  inner product
#  derivatives
#  laplacian
#  dotgrad

# zero(ej::FourierMode) = FourierMode(0, 0, ej.wavenumber)

FourierMode{T}() where {T<:Real} = FourierMode(zero(T),0,zero(T))
FourierMode() = FourierMode{T_ELEMENT}()
FourierMode(j, α) = FourierMode(one(α),j,α)

function (f::FourierMode)(x::Real) 
    if f.waveindex < 0 
        return f.coeff*cos(f.wavenumber*f.waveindex*x)
    elseif f.waveindex > 0
        return f.coeff*sin(f.wavenumber*f.waveindex*x)
    else
        return f.coeff
    end
end

compatible(ej::FourierMode, ek::FourierMode) = ej.wavenumber == ek.wavenumber

*(c::Real, ej::FourierMode) = FourierMode(c*ej.coeff, ej.waveindex, ej.wavenumber)

-(ej::FourierMode) = FourierMode(-ej.coeff, ej.waveindex, ej.wavenumber)

isorthogonal(ej::FourierMode, ek::FourierMode) = ej.waveindex != ek.waveindex 

# 1/Lx integral_0^Lx ej(x) ek(x) dx == 0 for j!=k, 1 for j=k=0, 1/2 for j=k!=0 (times coeffs)
function innerproduct(ej::FourierMode{T}, ek::FourierMode{T}) where {T<:Real}
    compatible(ej,ek) || error("incompatible ej=$(ej), ek=$(ek)")
    
    if isorthogonal(ej,ek) || ej.coeff == zero(T) || ek.coeff == zero(T) 
        return zero(T)

    elseif ej.waveindex == 0             # ek.waveindex == 0 as well, since ej, ek are nonorthogonal
        return ej.coeff*ek.coeff         # 1/Lx integral_0^Lx 1 dx = 1

    else
        return one(T)/2*ej.coeff*ek.coeff    # 1/Lx integral_0^Lx cos^2(ajx) dx = 1/Lx integral_0^Lx sin^2(ajx) dx = 1/2
    end
end

norm2(ej::FourierMode) = innerproduct(ej,ej)
norm(ej::FourierMode) = sqrt(norm2(ej))
    
# d/dx Ej(x) = c * (alpha*j) * E(-j)(x)
function derivative(ej::FourierMode)
    FourierMode(ej.coeff*ej.wavenumber*ej.waveindex, -ej.waveindex, ej.wavenumber)
end

# d^n/dx^n :
# d0/dx0  Ej(x) =                Ej(x) 
# d1/dx1  Ej(x) =  (alpha*j)   * E(-j)(x) 
# d2/dx2  Ej(x) = -(alpha*j)^2 * Ej(x)
# d3/dx3  Ej(x) = -(alpha*j)^3 * E(-j)(x)
# d4/dx4  Ej(x) =  (alpha*j)^4 * Ej(x),   etc
function derivative(ej::FourierMode, n::Int)
    n >= 0 || error("can't do nth derivative for negative n, n=$n")
    if n==0 
        return ej
    elseif n==1
        return FourierMode(ej.coeff*ej.wavenumber*ej.waveindex, -ej.waveindex, ej.wavenumber)
    elseif n==2
        return FourierMode(-ej.coeff*square(ej.wavenumber*ej.waveindex), ej.waveindex, ej.wavenumber)        
    else
        signj  = (n%2 == 0) ? 1 : -1 # implements above ±j pattern for E(±j)(x)
        signaj = (n%4 <  2) ? 1 : -1 # implements above ±j pattern for ±(alpha j)^n
        return FourierMode(ej.coeff*signaj*pow(ej.wavenumber*ej.waveindex, n), signj*ej.waveindex, ej.wavenumber)
    end
end

# multiplication of ej*ek = {cos(ajx) | 1 | sin(ajx)} * {cos(akx) | 1 | sin(akx)}
# via trig identities, e.g. cos(ajx)*cos(akx) = 1/2 cos(a(j-k)x) + 1/2 cos(a(j+k)x)
# return an array of FourierModes since some products have 1 term, some have 2. 
function *(ej::FourierMode{T}, ek::FourierMode{T}) where {T<:Real}
    !compatible(ej,ek) && error("incompatible ej, ek")

    # convenience precalculations, not used in all cases 
    # some could be pushed to just-in-time at the cost of clarity
    j = ej.waveindex
    k = ek.waveindex
    a = ej.wavenumber
    c = ej.coeff*ek.coeff
    c2 = T(c) / 2

    jksum =  abs(j) + abs(k)           # |j| + |k|
    jkdiff = abs(abs(j) - abs(k))      # ||j| - |k||
    jkdiffsign = sign(abs(j) - abs(k)) # sign(|j|-|k|)

    # The algorithm and indices here are so nonobvious!
    # The coding of FourierModes is 
    #   j<0 represents cos(a|j|x)
    #   j=0 represents 1
    #   j>0 represents sin(a|j|x)
    # So here's what the cases below represent. I'll leave out the coefficients cj, ck for simplicity
    # and leave absolute values in places where they could be simplified, for clarity of meaning.

    # case j==0: This represents ej ek == 1 x {cos(a|k|x), 1, or sin(a|k|x)}.
    # So the product ej ek == ek is given by ek with a modified coefficient.

    # case k==0: This represents ej ek == {cos(a|j|x), 1, or sin(a|j|x)} x 1.
    # So the product ej ek == ej is given by ej with a modified coefficient.

    # case j,k < 0: This represents ej ek == cos(a|j|x) cos(a|k|x). The product 
    # cos(a|j|x) cos(a|k|x) == 1/2 cos(a(|j|+|k|)x) + 1/2 cos(a(|j|-|k|)x) is then 
    # represented by a sum of FourierModes with -(|j|+|k|) == -jksum and -||j|-|k|| = -jkdiff

    # case j,k > 0. This represents ej ek == sin(a|j|x) sin(a|k|x). The product 
    # sin(a|j|x) sin(a|k|x) == 1/2 cos(a(|j|+|k|)x) - 1/2 cos(a(|j|-|k|)x) is then
    # represented by a sum of FourierModes with -(|j|+|k|) == -jksum and -||j|-|k|| = -jkdiff

    # case j + k == 0 (j,k != 0). We have j,k equal in magnitude with opposite sign,
    # representing ej ek == cos(a|j|x) sin(a|j|x). The product 
    # cos(a|j|x) sin(a|j|x) == 1/2 sin(a|2j|x)
    # is represented by a single FourierMode with waveindex -2|j|. 

    # case j<0 and k>0. This represents ej ek == cos(a|j|x) sin(a|k|x). The product
    # cos(a|j|x) sin(a|k|x) == 1/2 sin(a(|j|+|k|)x) - 1/2 sin(a(|j|-|k|)x)
    #                       == 1/2 sin(a(|j|+|k|)x) - 1/2 sin(a(||j|-|k||)x) sign(|j|-|k|)
    # is represented by two FourierModes with |j|+|k| = jksum and ||j|-|k|| = jkdiff.

    # case j>0 and k<0. This represents ej ek == sin(a|j|x) cos(a|k|x). The product
    # sin(a|j|x) cos(a|k|x) == 1/2 sin(a(|j|+|k|)x) + 1/2 sin(a(|j|-|k|)x) 
    #                       == 1/2 sin(a(|j|+|k|)x) + 1/2 sin(a(||j|-|k||)x) sign(|j|-|k|)
    # is represented by two FourierModes with |j|+|k| = jksum and ||j|-|k|| = jkdiff.


    if j == 0          # ej ek == 1 x ek == 1 x {cos(a|k|x), 1, or sin(a|k|x)}
        return [FourierMode(c, k, a)]

    elseif k == 0      # ej ek == ej x 1 == {cos(a|j|x), 1, or sin(a|j|x)} x 1
        return [FourierMode(c, j, a)]

    elseif j<0 && k<0  # ej ek == cos(a|j|x) cos(a|k|x) == 1/2 cos(a(|j|+|k|)x) + 1/2 cos(a(|j|-|k|)x)
        return [FourierMode(c2, -jkdiff, a); FourierMode(c2, -jksum, a);]
        
    elseif j>0 && k>0  # ej ek == sin(a|j|x) sin(a|k|x) == 1/2 cos(a(|j|+|k|)x) - 1/2 cos(a(|j|-|k|)x)
        return [FourierMode(c2,  -jkdiff, a); FourierMode(-c2, -jksum, a)]   
         
    elseif j + k == 0  # ej ek == cos(a|j|x) sin(a|j|x) == 1/2 sin(a|2j|x) (note j,k are nonzero)
        return [FourierMode(c2, abs(2j), a)]            

    elseif j<0 && k>0  # ej ek == cos(a|j|x) sin(a|k|x) == 1/2 sin(a(|j|+|k|)x) - 1/2 sin(a(|j|-|k|)x)
        return [FourierMode(c2, jksum, a); FourierMode(-jkdiffsign*c2, jkdiff, a)]            
        
    elseif j>0 && k<0  # ej ek == sin(a|j|x) cos(a|k|x) == 1/2 sin(a(|j|+|k|)x) + 1/2 sin(a(|j|-|k|)x)
        return [FourierMode(c2, jksum, a); FourierMode(jkdiffsign*c2, jkdiff, a)]            

    end

    error("oops, logic error, shouldn't get here in execution!")
    return [FourierMode(T(0), 0, T(0))]

end

############################################################################
"""
    BasisComponent{T<:Real}

BasisComponent represents a function of (x,y,z) that is Fourier in x,z and polynomial in y.
E.g. f(x,y,z) = c * E_j(x) * E_l(z) * P(y) for a real coefficient c, two Fourier modes E_j(x), E_l(z),
and a wall-normal polynomial p(y)
"""
struct BasisComponent{T<:Real}
    coeff::T
    ejx::FourierMode{T}    # represents cos(ajx), 1, or sin(ajx)
    ekz::FourierMode{T}    # represents cos(glz), 1, or sin(glz)
    p::Polynomial{T, :y}  # polynomial p(y) for wall-normal variation (rational coeffs)
    pparity::Int64      # parity of polynomial P(y): p==1 even, p==-1 odd, p==0 neither
    function BasisComponent(
        coeff::TC,
        ejx::FourierMode{TX},
        ekz::FourierMode{TZ},  
        p::Polynomial{TP, :y}, 
        pparity::Integer
    ) where {TC<:Real, TX<:Real, TZ<:Real, TP<:Real}
        T = promote_type(TC, TX, TZ, TP)
        return new{T}(
            T(coeff),
            FourierMode(
                T(ejx.coeff),
                ejx.waveindex,
                T(ejx.wavenumber)
            ),
            FourierMode(
                T(ekz.coeff),
                ekz.waveindex,
                T(ekz.wavenumber)
            ),
            Polynomial{T, :y}(p),
            Int64(pparity)
        )
    end
end

BasisComponent{T}() where {T<:Real} = BasisComponent(zero(T), FourierMode{T}(), FourierMode{T}(), Polynomial(zero(T), :y), 0) 
BasisComponent() = BasisComponent{T_ELEMENT}()

(f::BasisComponent)(x::Real, y::Real, z::Real) = f.coeff * f.ejx(x) * f.ekz(z) * f.p(y)
(f::BasisComponent)(x::Vector) = f.coeff * f.ejx(x[1]) * f.ekz(x[3]) * f.p(x[2])

zero(f::BasisComponent{T}) where {T<:Real} = zero(T)
zero(::Type{BasisComponent{T}}) where {T<:Real} = zero(T)

compatible(f::BasisComponent, g::BasisComponent) = compatible(f.ejx, g.ejx) && compatible(f.ekz, g.ekz)

regularize(f::BasisComponent) = BasisComponent(1, f.ejx, f.ekz, f.coeff*f.p)

# ors (||) between orthog checks because these terms are multiplicative
isorthogonal(f::BasisComponent, g::BasisComponent) = isorthogonal(f.ejx, g.ejx) || isorthogonal(f.ekz, g.ekz) || f.pparity*g.pparity  == -1 

# f = f.c f.ej(x) f.p(y) f.ek(z), g = g.c g.ej(x) g.p(y) g.ek(z)

# (f,g) = 1/(Lx Ly Lz) int_0^Lx int_0^Ly int_0^Lz  f g dx dy dx

# (f,g) = (f.c*g.c) (1/Lx f.ej(x) g.ej(x) dx) * (1/Ly int_-1^1 f.p(y) g.p(y) dy) * (1/Lz int_0^Lz f.ek(z) g.ek(z) dz)

# (f,g) = (f.c*g.c) (f.ej(x), g.ej(x)) * (1/Ly int_-1^1 f.p(y) g.p(y) dy) * (f.ek(z), g.ek(z))

function innerproduct(f::Polynomial{T}, g::Polynomial{T}, fparity=0, gparity=0) where {T<:Real}
    if fparity*gparity == -1
        return zero(T)
    end

    fgintegral = integrate(f*g)
    return T(fgintegral(1) - fgintegral(-1)) / 2
end

norm2(f::Polynomial) = innerproduct(f,f)
norm(f::Polynomial) = sqrt(norm2(f))

function innerproduct(f::BasisComponent{T}, g::BasisComponent{T}) where {T<:Real}
    compatible(f,g) || error("incompatible f,g")

    if isorthogonal(f,g) || f.coeff == 0 || g.coeff == 0 || f.pparity*g.pparity == -1
        return zero(T)
    end
    
    return f.coeff*g.coeff * innerproduct(f.ejx, g.ejx) * innerproduct(f.ekz,g.ekz) * innerproduct(f.p, g.p) 
end
    
norm2(f::BasisComponent) =  f.coeff^2 * norm2(f.ejx) * norm2(f.ekz) * norm2(f.p)
norm(f::BasisComponent) =  sqrt(norm2(f))
    
# df/dx_i, differentiate w.r.t ith coordinate of x, 
xderivative(f::BasisComponent) = BasisComponent(f.coeff, derivative(f.ejx), f.ekz, f.p,  f.pparity)  # d/dx
yderivative(f::BasisComponent) = BasisComponent(f.coeff, f.ejx, f.ekz, derivative(f.p), -f.pparity) # d/dy
zderivative(f::BasisComponent) = BasisComponent(f.coeff, f.ejx, derivative(f.ekz), f.p,  f.pparity)  # d/dx

xderivative(f::BasisComponent, n::Int) = BasisComponent(f.coeff, derivative(f.ejx,n), f.ekz, f.p, f.pparity)  # d/dx
yderivative(f::BasisComponent, n::Int) = BasisComponent(f.coeff, f.ejx, f.ekz, derivative(f.p,n), -f.pparity) # d/dy
zderivative(f::BasisComponent, n::Int) = BasisComponent(f.coeff, f.ejx, derivative(f.ekz,n), f.p, f.pparity)  # d/dx

*(c::Real, f::BasisComponent) = BasisComponent(c*f.coeff, f.ejx, f.ekz, f.p, f.pparity)
-(f::BasisComponent) = BasisComponent(-f.coeff, f.ejx, f.ekz, f.p, f.pparity)


function laplacian(f::BasisComponent)
    pyy = derivative(f.p, 2)
    aj2_gl2_p = (square(f.ejx.wavenumber*f.ejx.waveindex) + square(f.ekz.wavenumber*f.ekz.waveindex))*f.p
    BasisComponent(f.coeff, f.ejx, f.ekz, pyy - aj2_gl2_p, f.pparity)
end

function *(q::Polynomial, f::BasisComponent)
    pqparity = pparity*polyparity(q)
    BasisComponent(f.coeff, f.ejx, f.ekz, q*p, pqparity)
end

# return value is an array of BasisComponents, each elem representing a term in a sum
function *(f::BasisComponent, g::BasisComponent)

    fgcoeff = f.coeff*g.coeff

    # The products are arrays of FourierModes (with 1 to 4 elements), representing the terms of the product 
    fgejx = f.ejx * g.ejx  # {1 or cos(j1 x) or sin(j1 x)} * {1 or cos(j2 x) or sin(j2 x)} = sum of some trig functions
    fgekz = f.ekz * g.ekz  # {1 or cos(k1 z) or sin(k1 z)} * {1 or cos(k2 z) or sin(k2 z)} = sum of some trig functions

    #@show Estr.(fgejx,1)
    #@show Estr.(fgekz,2)

    # This product is a single polynomial 
    fgpoly = f.p * g.p
    fgpparity = f.pparity*g.pparity
    
    M = length(fgejx)
    N = length(fgekz)

    # Allocate and fill an array of BasisComponents to store up to four terms of form c ej(x) ek(z) p(y)
    rtn = fill(BasisComponent(), M*N)

    k = 1
    for m=1:M, n=1:N
        rtn[k] = BasisComponent(fgcoeff, fgejx[m], fgekz[n], fgpoly, fgpparity)
        k += 1
    end
    rtn
end

function innerproduct(f::BasisComponent{T}, g::AbstractVector{BasisComponent{T}}) where {T<:Real}
    ip = zero(T)
    for gi in g
        ip += innerproduct(f, gi)
    end
    return ip
        
end

function (f::Vector{BasisComponent{T}})(x::Real, y::Real, z::Real) where {T<:Real}
    rtn = zero(T)
    for fi in f
        rtn += fi(x, y, z)
    end
    return rtn
end


############################################################################

struct BasisFunction{T<:Real}
    u::SVector{3,BasisComponent{T}}
end

BasisFunction{T}() where {T<:Real} = BasisFunction(SVector(BasisComponent{T}(), BasisComponent{T}(), BasisComponent{T}()))
BasisFunction() = BasisFunction{T_ELEMENT}()

BasisFunction(u::BasisComponent, v::BasisComponent, w::BasisComponent) = BasisFunction(SVector(u, v, w))

(f::BasisFunction)(x::Real, y::Real, z::Real) = [f.u[1](x,y,z); f.u[2](x,y,z); f.u[3](x,y,z)]
(f::BasisFunction)(x::Vector) = [f.u[1](x); f.u[2](x); f.u[3](x)]

function compatible(f::BasisFunction, g::BasisFunction)
    compatible(f.u[1], g.u[1]) &&  compatible(f.u[2], g.u[2]) && compatible(f.u[3], g.u[3])
end

# ands (&&) between orthog checks because these terms sum in innerproduct
function isorthogonal(f::BasisFunction, g::BasisFunction)
    isorthogonal(f.u[1], g.u[1]) && isorthogonal(f.u[2], g.u[2]) && isorthogonal(f.u[3], g.u[3])
end

function innerproduct(f::BasisFunction, g::BasisFunction)
    # no use checking for f,g orthogonality, those same checks are done in componentwise innerproduct evals 
    ip1 = innerproduct(f.u[1], g.u[1])
    ip2 = innerproduct(f.u[2], g.u[2])
    ip3 = innerproduct(f.u[3], g.u[3])
    ip1 + ip2 + ip3
end

norm2(f::BasisFunction) = innerproduct(f,f)
norm(f::BasisFunction) = sqrt(norm2(f))

zero(psi::BasisFunction) = BasisFunction(zero.(psi.u))
xderivative(psi::BasisFunction, n::Int=1) = BasisFunction(xderivative(psi.u[1],n), xderivative(psi.u[2],n), xderivative(psi.u[3],n))
yderivative(psi::BasisFunction, n::Int=1) = BasisFunction(yderivative(psi.u[1],n), yderivative(psi.u[2],n), yderivative(psi.u[3],n))
zderivative(psi::BasisFunction, n::Int=1) = BasisFunction(zderivative(psi.u[1],n), zderivative(psi.u[2],n), zderivative(psi.u[3],n))

laplacian(psi::BasisFunction) = BasisFunction(laplacian.(psi.u))

*(c::Real, f::BasisFunction) = BasisFunction(c*f.u)
-(f::BasisFunction) = BasisFunction(-1*f.u)

function *(q::Polynomial, f::BasisFunction)
    qparity = polyparity(q)
    u = BasisComponent(f.u[1].coeff, f.u[1].ejx, f.u[1].ekz, q*f.u[1].p, qparity*f.u[1].pparity)
    v = BasisComponent(f.u[2].coeff, f.u[2].ejx, f.u[2].ekz, q*f.u[2].p, qparity*f.u[2].pparity)
    w = BasisComponent(f.u[3].coeff, f.u[3].ejx, f.u[3].ekz, q*f.u[3].p, qparity*f.u[3].pparity)
    BasisFunction(u,v,w)
end
    
"""
For psi = [u,v,w], return [v,0,0]
"""
function vex(psi::BasisFunction{T}) where {T<:Real}
    u = BasisComponent(psi.u[2].coeff, psi.u[2].ejx, psi.u[2].ekz, psi.u[2].p, psi.u[2].pparity)
    v = BasisComponent(zero(T), psi.u[2].ejx, psi.u[2].ekz, Polynomial(zero(T), :y), 0)
    w = BasisComponent(zero(T), psi.u[3].ejx, psi.u[3].ekz, Polynomial(zero(T), :y), 0)
    BasisFunction(u,v,w)
end


# returns a 3-vector of vectors of BasisComponents
# fdotgradg[1] = vector of BasisComponents representing sum ejx ekz p(y) for u component
# fdotgradg[2] = vector of BasisComponents representing sum ejx ekz p(y) for v component
# fdotgradg[3] = vector of BasisComponents representing sum ejx ekz p(y) for w component

function dotgrad(f::BasisFunction{T}, g::BasisFunction{T}) where {T<:Real}
    rtn = Vector{Vector{BasisComponent{T}}}(undef, 3)    
    rtn[1] = [f.u[1]*xderivative(g.u[1]); f.u[2]*yderivative(g.u[1]); f.u[3]*zderivative(g.u[1])]
    rtn[2] = [f.u[1]*xderivative(g.u[2]); f.u[2]*yderivative(g.u[2]); f.u[3]*zderivative(g.u[2])]
    rtn[3] = [f.u[1]*xderivative(g.u[3]); f.u[2]*yderivative(g.u[3]); f.u[3]*zderivative(g.u[3])]

    rtn
end

function innerproduct(
    f::BasisFunction{T}, 
    g::AbstractVector{V}) where {T<:Real, V<:AbstractVector{<:BasisComponent{T}}}
    s = zero(T)
    for i=1:3  # for each component f[i], g[i]
        for n=1:length(g[i])
            s += innerproduct(f.u[i], g[i][n]) # sum the inner products of f[i] with terms in g[i]
        end
    end
    s
end

############################################################################
# Utility functions from here on

"""
make an array of indices -J <= j <= J in order [0, 1, -1, 2, -2, ..., J, -J]
"""
function fourierIndices(J)
    jset = fill(0, 2J+1)
    jset[2:2:end] = -(1:J)
    jset[3:2:end] = 1:J
    jset
end

"""
    basisIndices(J::Int, K::Int, L::Int)
    basisIndices(J::Int, K::Int, L::Int, H::Vector{Symmetries})

Return the set of i,j,k,l basis set indices for J,K,L discretization limits with -J ≤ j ≤ J,
-K ≤ k ≤ K, and 0 ≤ l ≤ L, and 1 ≤ i ≤ 6. The allowed values of j,k,l depend on i. Return
value is an N x 4 matrix each row of which specifies an i,j,k,l index for a basis element.

If symmetry subgroup H is specified, filter the set of indices to include only those that
are symmetric in all the symmetries in H. (H should be a vector of generators for the group).
"""
function basisIndices(J::Int, K::Int, L::Int)
    N = (2J+1)*(2K+1)*(2L+1) + 1
    ijkl = fill(0, (N,4))
    Jset = fourierIndices(J)
    Kset = fourierIndices(K)
    n = 1
    for j in Jset, k in Kset

        # assign all \Psi_ijkl for the i and l values reachable by this j,k pair
        
        if j == 0
            # insert L+1 type i==1 elements Psi_1jkl = [E_k(z) S_l'(y); 0; 0] with given k and l = 0:L,
            ijkl[n:n+L, 1] .= 1   # assign i==1 
            ijkl[n:n+L, 2] .= 0   # assign j==0
            ijkl[n:n+L, 3] .= k   # assign k
            ijkl[n:n+L, 4] .= 0:L # assign l == 0:L
            n += L+1
        end

        if j == 0 && k != 0
            # insert L type i==4 elements Psi_2jkl = [0; γk E_k(z) S_l(y); E_{-k}(z) S_l'(y)] with given k and l = 1:L,        
            ijkl[n:n+L-1, 1] .= 2   # assign i==2 
            ijkl[n:n+L-1, 2] .= 0   # assign j==0
            ijkl[n:n+L-1, 3] .= k   # assign k
            ijkl[n:n+L-1, 4] .= 1:L # assign l == 1:L
            n += L
        end
        
        if k == 0   
            # insert L+1 type i==3 elements Psi_3jkl = [0; 0; E_j(x) R_l(y)] with the given j and l = 0:L
            ijkl[n:n+L, 1] .= 3   # assign i==3
            ijkl[n:n+L, 2] .= j   # assign j
            ijkl[n:n+L, 3] .= 0   # assign k==0
            ijkl[n:n+L, 4] .= 0:L # assign l == 0:L
            n += L+1
        end

        if k == 0 && j != 0
            # insert L type i==4 elements Psi_4jkl = [E_{-j}(x) S_l'(y); αj E_j(x) S_l(y); 0] with the given j and l = 1:L
            ijkl[n:n+L-1, 1] .= 4   # assign i==4 
            ijkl[n:n+L-1, 2] .= j   # assign j
            ijkl[n:n+L-1, 3] .= 0   # assign k==0
            ijkl[n:n+L-1, 4] .= 1:L # assign l == 0:L
            n += L
        end
        
        # if j and k are nonzero, then there will be
        # L+1 type 5 basis elements Psi_5jkl with the given j,k values and l = 0:L
        # L   type 6 basis elements Psi_6jkl with the given j,k values and l = 1:L
        #
        # Psi_5jkl = [ γk E_{-j}(x) E_k(z) S_l'(y);
        #              0;
        #             -αj E_j(x) E_{-k}(z) S_l'(y) ]
        #
        # Psi_6jkl = [ γk   E_{-j}(x) E_k(z) S_l'(y);
        #              αγjk E_j(x)    E_k(z) S_l(y); 
        #              αj   E_j(x) E_{-k}(z) S_l'(y) ]
        if j != 0  && k != 0

            # insert L+1 type 5 basis elements Psi_5jkl with the given j,k values and l = 0:L
            ijkl[n:n+L, 1] .= 5    # assign i==5
            ijkl[n:n+L, 2] .= j    # assign j
            ijkl[n:n+L, 3] .= k    # assign k
            ijkl[n:n+L, 4] .= 0:L  # assign l == 0:L
            n += L+1

            # insert L type 6 basis elements Psi_6jkl with the given j,k values and l = 1:L
            ijkl[n:n+L-1, 1] .= 6     # assign i==4
            ijkl[n:n+L-1, 2] .= j     # assign j==j
            ijkl[n:n+L-1, 3] .= k     # assign k==j
            ijkl[n:n+L-1, 4] .= 1:L   # assign l == 1:L
            n += L
        end
    end
    ijkl
end

function basisIndices(J::Int, K::Int, L::Int, H::AbstractVector{Symmetry})
    # make full set of indices with no symmetry restrictions
    ijklfull = basisIndices(J,K,L)

    # get indices of all basis elements that are symmetric in all generators
    nsymm = findall([symmetric(ijklfull[n,:], H) for n in 1:size(ijklfull,1)])

    ijklfull[nsymm, :]
end


"""
return OffSetArray of Legendre polynomials P_n(y) for n=0 to K
"""
function legendrePolynomials(T::Type, L::Int, symbol=:y)
    L >= 0 || error("invalid L value L=$(L), should be nonnegative")
    P = OffsetVector(fill(Polynomial([zero(T)], :y), L+1), -1)
    P[0] = Polynomial{T}([1], :y)
    if L>1
        P[1] = Polynomial{T}([0, 1], :y)
    end
    for n = 2:L
        P[n] = ((2*n-1)*P[1]*P[n-1] - (n-1)*P[n-2])/n
    end
    return P
end
legendrePolynomials(L, symbol=:y) = legendrePolynomials(T_ELEMENT, L, symbol)


"""
    ratsqrt(n::Int)

Return a rational approximation to sqrt(n). Useful for normalizing 
rational-coeff polynomials while staying rational, not converting to float.
"""
function ratsqrt(n::Int)
    r = Int(round(sqrt(n)))
    if r^2 == n
        return r
    end
    r = 1//2 * (r + n//r)
end

"""
    Snormalizer(n)

For functions Sn(y) = (1-y^2)^2 P_{n-1}(y), n >= 1, and S0(y) = (1-y^2)y,
norm(Sn) ~ 1/sqrt(n) and norm(Sn') ~ sqrt(n). Snormalizer is a rational-valued
function that scales as Snormalizer(n) ~ 1/sqrt(n), so that 
norm(Sn/Snormalizer(n)) ~ 1 and norm(Sn'/Snormalizer(n)) ~ n.
"""
function Snormalizer(n)
    n == 0 && return 1//2
    return 1//2 / ratsqrt(n)
end

min1(n) = n < 1 ? 1 : n

"""
make a set of basis functions Ψijkl, return as an array of BasisFunctions
   α,γ == Fourier wavenumbers 2π/Lx, 2π/Lz, 
   basis functions are all zero divergence, div Psi = 0, and zero at walls, Psi(x,±1,z) = 0.
   set definition is complicated, see LaTeX documentation elsewhere...
"""
function basisSet(α::T, γ::T, ijkl::Matrix{Int}; normalize=false) where {T<:Real}
    N = size(ijkl, 1)
    J = maximum(abs.(ijkl[:,2]))
    K = maximum(abs.(ijkl[:,3]))
    L = maximum(ijkl[:,4])

    smasher = Polynomial{T}([1; 0; -1], :y)  # smasher  =  1-y^2
    smasher2 = smasher*smasher            # smasher2 = (1-y^2)^2
    P = legendrePolynomials(T, L)            # l = 0:L, even/odd with l

    S0 = Polynomial{T}([0; 1; 0; -1/3], :y)           # S0(y) = y - y^3/3
    S = OffsetVector([S0; smasher2.*P[0:L-1]], -1) # l = 1:L, even/odd with l+1
    Sprime = derivative.(S)                        # l = 0:L, even/odd with l

    Ejx = OffsetVector([FourierMode(j, α) for j in -J:J], -(J+1))
    Ekz = OffsetVector([FourierMode(k, γ) for k in -K:K], -(K+1))

    zeroEjx = FourierMode(zero(T), 0, α)
    zeroEkz = FourierMode(zero(T), 0, γ)
    zeropoly = Polynomial(zero(T), :y)
    zerocomp = BasisComponent(0, zeroEjx, zeroEkz, zeropoly, 0)

    # Note: the normalizer constants commented out below are an attempt to do
    # approximate normalization of rational BasisFunctions.
    Ψ = fill(BasisFunction(), N)
    for n=1:N
        i,j,k,l = ijkl[n,:]
        normalizer = one(T)
        if i == 1       # j=0
            #normalizer = normalize ? T(min1(l)*Snormalizer(l)) : one(T)
            Ψu = BasisComponent(one(T), Ejx[0], Ekz[k], Sprime[l], parity(l))
            Ψ[n] = BasisFunction(Ψu, zerocomp, zerocomp)

        elseif i == 2   # j=0, k≠0, l≠0
            #normalizer = normalize ? T((γ*abs(k) + min1(l))*Snormalizer(l)) : one(T)
            Ψv = BasisComponent(γ*k, Ejx[0], Ekz[k], S[l], parity(l+1))
            Ψw = BasisComponent(one(T),   Ejx[0], Ekz[-k], Sprime[l], parity(l))
            Ψ[n] = BasisFunction(zerocomp, Ψv, Ψw)

        elseif i == 3   # k=0
            #normalizer = normalize ? T(min1(l)*Snormalizer(l)) : one(T)
            Ψw = BasisComponent(one(T), Ejx[j], Ekz[0], Sprime[l], parity(l))
            Ψ[n] = BasisFunction(zerocomp, zerocomp, Ψw)

        elseif i == 4   # j≠0, k=0, l≠0
            #normalizer = normalize ? T((l + α*abs(j))*Snormalizer(l)) : one(T)
            Ψu = BasisComponent(one(T), Ejx[-j], Ekz[0], Sprime[l], parity(l))
            Ψv = BasisComponent(α*j, Ejx[j], Ekz[0], S[l], parity(l+1))
            Ψ[n] = BasisFunction(Ψu, Ψv, zerocomp)

        elseif i == 5   # j≠0, k≠0
            #normalizer = normalize ? T((α*abs(j) + γ*abs(k))*min1(l)*Snormalizer(l)) : one(T)
            Ψu = BasisComponent( γ*k, Ejx[-j], Ekz[k], Sprime[l], parity(l)) 
            Ψw = BasisComponent(-α*j, Ejx[j], Ekz[-k], Sprime[l], parity(l))
            Ψ[n] = BasisFunction(Ψu, zerocomp, Ψw)

        elseif i == 6   # j≠0, k≠0, l≠0
            #normalizer = normalize ? T((l*γ*abs(k) + 2α*γ*abs(j*k) + l*α*abs(j))*Snormalizer(l)) : one(T)
            Ψu = BasisComponent(γ*k,      Ejx[-j], Ekz[k], Sprime[l], parity(l))
            Ψv = BasisComponent(2α*γ*j*k, Ejx[j], Ekz[k], S[l], parity(l+1))                
            Ψw = BasisComponent(α*j,      Ejx[j], Ekz[-k], Sprime[l], parity(l))
            Ψ[n] = BasisFunction(Ψu, Ψv, Ψw)

        else
            println("oops! logic error, execution shouldn't get here; i==$(i)")
            Ψ[n] = BasisFunction(zerocomp, zerocomp, zerocomp)
        end
        if (normalize)
            Ψ[n] = 1/norm(Ψ[n]) * Ψ[n]
        end
    end
    Ψ
end



"""
make a set of basis functions Ψijkl 
   α,γ == Fourier wavenumbers 2π/Lx, 2π/Lz, 
   J,K,L == discretization limits: -J ≤ j ≤ J, -K ≤ k ≤ K, 0 ≤ l ≤ L.
            (j,k, are x,z Fourier indices, l is y polynomial index, i is 1 through 6 for diff forms of Ψ
   basis functions are all zero divergence, div Psi = 0, and zero at walls, Psi(x,±1,z) = 0.
"""
function basisSet(α::Real, γ::Real, J::Int, K::Int, L::Int; normalize=false)
    ijkl = basisIndices(J,K,L)
    Ψ = basisSet(α, γ, ijkl, normalize=normalize)
end

m2xyz(m) = m==1 ? "x" : (m==2 ? "y" : "z") # produce string x,y, or z from integer index m 
m2αγ(m)  = m==1 ? "α" : "γ"                # produce string α or γ from integer index m 

function coeffstring(x) 
    if x ≈ 1 
        return "(+1)"
    elseif x ≈ -1
        return "(-1)"
    else
        return "($(x))"
    end
end

# produce string "sin(α x)", "cos(γ z)"" etc from integer indices
function Estr(j,m) 
    if j < -1
        return "ACK cos($(abs(j))$(m2αγ(m))$(m2xyz(m)))"
    elseif j == -1
        return "ACK cos($(m2αγ(m))$(m2xyz(m)))"
    elseif j==0
        return "ACK 1"
    elseif j==1
        return "ACK sin($(m2αγ(m))$(m2xyz(m)))"
    else        
        return "ACK sin($(j)$(m2αγ(m))$(m2xyz(m)))"
    end
end

# produce string "2 sin(α x) cos(γ z)" etc 
function Estr(E::FourierMode, m::Int) 
    if E.coeff == 0
        return "0 "
    end

    cstr = coeffstring(E.coeff)
    j = E.waveindex

    if j < -1
        return "$(cstr) cos($(abs(j))$(m2αγ(m))$(m2xyz(m)))"
    elseif j == -1
        return "$(cstr) cos($(m2αγ(m))$(m2xyz(m)))"
    elseif j==0
        return "$(cstr)"
    elseif j==1
        return "$(cstr) sin($(m2αγ(m))$(m2xyz(m)))"
    else  # j>1
        return "$(cstr) sin($(j)$(m2αγ(m))$(m2xyz(m)))"
    end
end

function estr(E::FourierMode)
    if E.coeff == 0
        return "0 "
    end

    cstr = coeffstring(E.coeff)
    j = E.waveindex
    m = 1
    if j < -1
        return "$(cstr) cos($(abs(j))$(m2αγ(m))$(m2xyz(m)))"
    elseif j == -1
        return "$(cstr) cos($(m2αγ(m))$(m2xyz(m)))"
    elseif j==0
        return "$(cstr)"
    elseif j==1
        return "$(cstr) sin($(m2αγ(m))$(m2xyz(m)))"
    else  # j>1
        return "$(cstr) sin($(j)$(m2αγ(m))$(m2xyz(m)))"
    end
end

function ustr(u::BasisComponent)
    if u.coeff == 0
        return "0"
    else 
        return "$(coeffstring(u.coeff)) $(Estr(u.ejx, 1)) $(Estr(u.ekz, 3)) ($(u.p)))"
    end
end

function ustr(u::Vector{BasisComponent})
    rtn = "( "
    for i in 1:length(u)
        if u[i].coeff != 0
            rtn *= " + $(coeffstring(u[i].coeff)) $(Estr(u[i].ejx, 1)) $(Estr(u[i].ekz, 3)) ($(u[i].p))) \n"
        end
    end
    rtn *= " )"
end

#===========================================================================
function ustr(u::BasisComponentSum)
    
    sb = StringBuilder()
    append!(sb, "(");
    
    for i = 1:length(u.ejx)
        append!(sb, "($(Estr(u.ejx[i], 1)) $(Estr(u.ekz[i], 3)))")
        if i<length(u.ejx)
            append!(sb, " + ")
        end
    end

    append!(sb, ")");
    append!(sb, "($(u.p)))")
    
    return String(sb)
end

function psistr(psi::BasisFunctionSum)
    return "[  $(ustr(psi.u[1])) ; $(ustr(psi.u[2])) ; $(ustr(psi.u[3])) ]"

===========================================================================#

function psistr(psi::BasisFunction)
    return "[  $(ustr(psi.u[1])) ; $(ustr(psi.u[2])) ; $(ustr(psi.u[3])) ]"
end

# The logic of the following functions is best understood by comparing the if-else statements
# on the indices to the math expressions for Ψijkl(x,y,z), which have three components each built
# from
#   Ej(x) = cos ajx, 1, or sin ajx,
#   Ek(z) = cos gkz, 1, or sin gkz, and
#   Pl(y) = even/odd polynomials
# depending on the ijkl indices. The symmetries of Ψijkl(x,y,z) are determined by the choice
# of cos/sin/1 in x,z and polynomial parities.

"""
return ±1 for Ψijkl symmetric/antisymmetric under x reflection σx
"""
function xreflection(ijkl::Vector{Int})
    i,j,k,l = ijkl
    #   i values | rtn value
    #   1        | -1
    #   2        |  1
    #   3 4 5 6  | (j <= 0) ? 1 : -1
    if i == 1 
        return -1
    elseif i == 2
        return 1
    elseif  3 <= i <= 6
        return (j <= 0) ? 1 : -1
    else
        error("invalid index i = $(i) not in 1:6")
    end
end

"""
return ±1 for Ψijkl symmetric/antisymmetric under y reflection σy
"""
function yreflection(ijkl::Vector{Int})
    i,j,k,l = ijkl
    return l%2 == 0 ? 1 : -1
end

"""
return ±1 for Ψijkl symmetric/antisymmetric under z translation τz
"""
function zreflection(ijkl::Vector{Int})
    i,j,k,l = ijkl
    #   i values | rtn value    
    #   3        | -1
    #   4        |  1
    #   1 2 5 6  | (k <= 0) ? 1 : -1
    if i == 3 
        return -1
    elseif i == 4
        return 1
    elseif i == 1 || i == 2 || i == 5 || i == 6
        return (k <= 0) ? 1 : -1
    else
        error("invalid index i = $(i) not in 1:6")
    end
end

"""
return ±1 for Ψijkl symmetric/antisymmetric under x translation by Lx/2
"""
function xtranslationLx2(ijkl::Vector{Int})
    i,j,k,l = ijkl
    return j%2 == 0 ? 1 : -1
end


"""
return ±1 for Ψijkl symmetric/antisymmetric under z translation by Lz/2
"""
function ztranslationLz2(ijkl::Vector{Int})
    i,j,k,l = ijkl
    return k%2 == 0 ? 1 : -1
end   

"""
return 1 for if |j|<=J, k<=K, |l|<=L
"""
function loworder(ijkl::Vector{Int}, JKL::Vector{Int})
    i,j,k,l = ijkl
    J,K,L = JKL
    return abs(j) <= J && abs(k) <= K && 0 <= l <= L
end   

xreflection(ijkl::Matrix{Int}) =  [xreflection(ijkl[n,:]) for n=1:size(ijkl,1)]
yreflection(ijkl::Matrix{Int}) =  [yreflection(ijkl[n,:]) for n=1:size(ijkl,1)]
zreflection(ijkl::Matrix{Int}) =  [zreflection(ijkl[n,:]) for n=1:size(ijkl,1)]
xtranslationLx2(ijkl::Matrix{Int}) =  [xtranslationLx2(ijkl[n,:]) for n=1:size(ijkl,1)]
ztranslationLz2(ijkl::Matrix{Int}) =  [ztranslationLz2(ijkl[n,:]) for n=1:size(ijkl,1)]
loworder(ijkl::Matrix{Int}, JKL::Vector{Int}) = [loworder(ijkl[n,:], JKL) for n=1:size(ijkl,1)]

function ijkl2file(ijkl, filebase; comment_char='#')
    filename = occursin(".asc", filebase) ? filebase : filebase * ".asc"
    io = open(filename, "w")
    M, N = size(ijkl)
    N == 4 || error("ijkl matrix should have 4 cols, but it has N=$N")
    println(io, "$comment_char $M")
    for n=1:M
        #i,j,k,l = ijkl[n,:]
        println(io, "$(ijkl[n,1]) $(ijkl[n,2]) $(ijkl[n,3]) $(ijkl[n,4])")
    end
    close(io)
end

function save(A::Matrix, filebase; comment_char='#')
    filename = occursin(".asc", filebase) ? filebase : filebase * ".asc"
    io = open(filename, "w")
    M,N = size(A)
    println(io, "$comment_char $M $N")
    for i=1:M, j=1:N
        print(io, A[i,j], j<N ? ' ' : '\n')
    end
    close(io)
end

function save(x::Vector, filebase; comment_char='#')
    filename = occursin(".asc", filebase) ? filebase : filebase * ".asc"
    io = open(filename, "w")
    N = length(x)
    println(io, "$comment_char $N")
    for i=1:N
        print(io, x[i], '\n')
    end
    close(io)
end


#=======================================
"""
    shear(x) = shearFunction(Ψ)

return a function shear(x) that computes the wall shear rate of u = sum x[i] Ψ[i]
"""
function shearFunction(Ψ::AbstractVector{BasisFunction{T}}) where {T<:Real}
    # Produce vector shearΨ for evaluating dissipation from x values:
    N = length(Ψ)
    println("Building shear(x)...")
    shearΨ = zeros(T, N)
    for i=1:N    
        if Ψ[i].u[1].ejx.waveindex == 0 && Ψ[i].u[1].ekz.waveindex == 0
            dΨudy = yderivative(Ψ[i].u[1])
            shearΨ[i] = Ψ[i].u[1].coeff*(dΨudy.p(1.0)  + dΨudy.p(-1.0))/2
        end
    end
    shear(x) = dot(shearΨ, x) + 1.0
end
======================================#

"""
    basis_index_dict(ijklΨ, ijklΦ)

Return n->m dict such that ijklΨ[n] = ijklΦ[dict[n],:] for all keys n in dict
"""
function basis_index_dict(ijklΨ, ijklΦ)
    n2m = Dict{Int,Int}([])
    for n = 1:size(ijklΨ,1)
        m = findfirst(row -> row == ijklΨ[n,:], eachrow(ijklΦ))  # find m s.t. ijklΦ[m,:] == ijklΨ[n,:]
        if m != nothing
            n2m[n] = m
        end
    end
    n2m
end


"""
    changebasis(x, ijklΨ, ijklΦ)
    
Convert expansion coefficients x in basis Ψ to expansion coeffs y in basis Φ.
so sum_n x[n] Ψ[n] = sum_m y[m] Φ[m] as closely as possible. Do this by
identifying elements of Φ that occur in Ψ and copying coefficients to match. 
"""
function changebasis(x, ijklΨ, ijklΦ)
    NmodesΨ = size(ijklΨ,1)
    NmodesΦ = size(ijklΦ,1)
    # input is sum_n x[n] Ψ[n], want to reexpress as sum_n y[n] Φ[n]
    # where Φ is superset of Ψ. (If Φ is subset of Ψ, drop missing modes.)
    # n is index into basis Ψ, m is index into basis Φ
    # n2m maps n in m s.t. ijklΦ[n2m[n]] = ijklΨ[n] and y[n2m[n]] = x[n]
    n2m = basis_index_dict(ijklΨ, ijklΦ) 
    y = zeros(NmodesΦ)
    for n in 1:NmodesΨ
        if haskey(n2m, n)
            y[n2m[n]] = x[n]
        end
    end
    y
end


# ====================================================================================
# functions for generating symmetric basis sets and ODE models



