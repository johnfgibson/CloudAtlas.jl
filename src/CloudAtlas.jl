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

# Visualization functions - these are only available when CairoMakie is loaded
# The actual implementations are in ext/CloudAtlasVisualizationExt.jl

# Note: we define PlotSettings here, because extensions cannot nicely export strutures.
"""
    PlotSettings

Configuration for velocity field plots.
"""
Base.@kwdef struct PlotSettings
    num_points::Int = 30
    arrow_scale::Float64 = 0.5
    arrow_lengthscale::Float64 = 0.8
    arrow_tiplength::Float64 = 14.0   # Values in pixels 
    arrow_tipwidth::Float64 = 14.0    # Values in pixels
    arrow_shaftwidth::Float64 = 2.0   # Value in pixels
    colormap::Vector{Symbol} = [:navyblue, :aqua, :lime, :orange, :red4]
    fig_size::Tuple{Int,Int} = (900, 1200)
end

function velocity_fields end
function velocity_fields_dns end
function velocity_fields_comparison end
export velocity_fields, velocity_fields_dns, velocity_fields_comparison, PlotSettings

end # module CloudAtlas
