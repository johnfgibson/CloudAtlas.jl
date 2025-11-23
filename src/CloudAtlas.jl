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

include("Symmetries.jl")

export Symmetry, symmetric, halfbox_symmetries

include("SparseBilinear.jl")

export SparseBilinear, sparse 

include("BasisFunctions.jl")

export FourierMode, BasisComponent, BasisFunction, compatible, isorthogonal, innerproduct, derivative, xderivative, yderivative, zderivative, *, zero, regularize, laplacian, dotgrad, fourierIndices, basisIndices, basisSet, estr, Estr, ustr, psistr, legendrePolynomials, xreflection, yreflection, zreflection, xtranslationLx2, ztranslationLz2, vex, norm, norm2, loworder, ijkl2file, save, polyparity, basis_index_dict, changebasis

include("ODEModels.jl")

export ODEModel, shear, length

include("Hookstep.jl")

export hookstepsolve, SearchParams

include("TWModels.jl")

export TWModel, has_shift_symmetry, save_sigma

end # module CloudAtlas
