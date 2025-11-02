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

include("BasisFunctions.jl")

export FourierMode, BasisComponent, BasisFunction, compatible, isorthogonal, innerproduct, derivative, xderivative, yderivative, zderivative, *, zero, regularize, laplacian, dotgrad, fourierIndices, basisIndexMap, basisIndexMap2, makeBasisSet, estr, Estr, ustr, psistr, legendrePolynomials, xreflection, yreflection, zreflection, xtranslationLx2, ztranslationLz2, intersect, findsymmetric, vex, norm, norm2, loworder, ijkl2file, save, SparseBilinear, sparse, polyparity

include("Hookstep.jl")

export hookstepsolve 

# include("Symmetries.jl")
# include("ODEModels.jl")

end # module CloudAtlas
