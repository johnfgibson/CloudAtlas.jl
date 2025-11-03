module CloudAtlas

using Polynomials
using StaticArrays
using OffsetArrays
using LinearAlgebra
using SparseArrays
using StringBuilders
import Base: *, -, zero, intersect
import LinearAlgebra: norm
import Polynomials: derivative
import SparseArrays: sparse

include("Hookstep.jl")

export hookstepsolve 

include("Symmetries.jl")

export Symmetry, symmetric

include("BasisFunctions.jl")

export FourierMode, BasisComponent, BasisFunction, compatible, isorthogonal, innerproduct, derivative, xderivative, yderivative, zderivative, *, zero, regularize, laplacian, dotgrad, fourierIndices, basisIndices, basisSet, estr, Estr, ustr, psistr, legendrePolynomials, xreflection, yreflection, zreflection, xtranslationLx2, ztranslationLz2, vex, norm, norm2, loworder, ijkl2file, save, SparseBilinear, sparse, polyparity, shearFunction

include("ODEModels.jl")

export ODEModel

end # module CloudAtlas
