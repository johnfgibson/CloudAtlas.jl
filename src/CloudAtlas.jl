module CloudAtlas

export greet

"Write a friendly message"
greet() = print("hello, world! (revised)")

include("Hookstep.jl")
include("BasisFunctions.jl")
include("Symmetries.jl")
include("ODEModels.jl")

end # module CloudAtlas
