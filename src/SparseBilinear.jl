# Sparse representation of bilinear operator N_ijk for calculating y_i = N_ijk x_j z_k or N_ijk x_j x_k

struct SparseBilinear{T} 
    ijk::Matrix{Int64}
    val::Vector{T}
    m::Int   # output dimension
    function SparseBilinear(ijk, val::Vector{T}, m) where T <: Real
        size(ijk,2) == 3 || error("ijk must have three columns")
        size(ijk,1) == length(val) || error("ijk and val must have same number of rows")
        new{T}(ijk, val, m)
    end
    
end

"""
construct SparseBilinear operator from vectors I,J,K of indices and m-vector V of values
"""
SparseBilinear(I, J, K, V::Vector{T}, m) where T<:Number = SparseBilinear{T}([I J K], V, m)


"""
construct SparseBilinear from dense Array{T,3}
"""
function SparseBilinear(N::Array{T, 3}) where T<:Number 
    m = size(N,1)
    # count number nonzeros
    nnz = 0
    for i=1:size(N,1), j=1:size(N,2), k=1:size(N,3)
        if N[i,j,k] != 0
            nnz += 1
        end
    end

    ijk = zeros(Int, nnz, 3)
    val = zeros(T, nnz)

    r = 1
    for i=1:size(N,1), j=1:size(N,2), k=1:size(N,3)
        if N[i,j,k] != 0
            ijk[r,1] = i
            ijk[r,2] = j
            ijk[r,3] = k
            val[r] = N[i,j,k]
            r += 1
        end
    end
    
    SparseBilinear(ijk, val, m)
end

sparse(N::Array{T, 3}) where T<:Number = SparseBilinear(N)

"""
evaluate N_ijk x_j y_k
"""
function (N::SparseBilinear{TN})(x::AbstractVector{TXY}, y::AbstractVector{TXY}) where {TN<:Real, TXY<:Real}
    T = promote_type(TN, TXY)
    rtn = zeros(T, N.m)
    for r = 1:length(N.val)
        # clearest form:  i,j,k = N.ijk[r,:]; rtn[i] += N.val[r] * x[j] * y[k]
        # same without temporaries:
        rtn[N.ijk[r,1]] += N.val[r] * x[N.ijk[r,2]] * y[N.ijk[r,3]]
    end
    rtn
end

"""
evaluate N_ijk x_j y_k
"""
(N::SparseBilinear{T})(x) where T<:Number = N(x,x)

"""
evaluate DN_ij = d/dx_j (N_ilk x_l x_k) 
"""
function derivative(N::SparseBilinear{TN}, x::AbstractVector{TX}) where {TN<:Real, TX<:Real}
    T = promote_type(TN, TX)
    DN = zeros(T, N.m, N.m) # use a dense matrix, since DN is 99% dense
    
    for r = 1:length(N.val)
        #println("$r "); flush(stdout)
        # clearest form
        # i,j,k = ijk[r,:]
        # DN[i, j] += N.val[r]*x[k]
        # DN[i, k] += N.val[r]*x[j]

        # same without temporaries
        DN[N.ijk[r,1], N.ijk[r,2]] += N.val[r] * x[N.ijk[r,3]]
        DN[N.ijk[r,1], N.ijk[r,3]] += N.val[r] * x[N.ijk[r,2]]
    end
    DN
end

function save(A::SparseMatrixCSC, filebase)
    filename = occursin(".asc", filebase) ? filebase : filebase * ".asc"
    io = open(filename, "w")
    M,N = size(A)
    Nnonzero = nnz(A)
    println(io, "% $Nnonzero 3")
    for j in 1:N
        for n in A.colptr[j]:(A.colptr[j+1]-1) # location in rowval, nzval where col j data is stored
            println(io, "$(A.rowval[n]) $j $(A.nzval[n])")
        end
    end
    close(io)
end

