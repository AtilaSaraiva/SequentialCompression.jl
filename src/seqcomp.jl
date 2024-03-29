"""
    CompressedArraySeq{T,Nx}

A mutable structure for storing time-dependent arrays in a compressed format.

# Fields
- `data::Vector{UInt8}`: Compressed data in byte form.
- `headpositions::Vector{Int64}`: Positions of the beginning of each time slice in `data`.
- `tailpositions::Vector{Int64}`: Positions of the end of each time slice in `data`.
- `spacedim::NTuple{Nx,Int32}`: Dimensions of the spatial grid.
- `timedim::Int32`: Number of time steps.
- `eltype::Type{T}`: Element type of the uncompressed array.
- `tol::Float32`: [Mean absolute error that is tolerated](https://zfp.readthedocs.io/en/release0.5.5/modes.html#fixed-accuracy-mode).
- `precision::Float32`: [Controls the precision, bounding a weak relative error](https://zfp.readthedocs.io/en/release0.5.5/modes.html#fixed-precision-mode).
- `rate::Int64`: [Fixes the bits used per value](https://zfp.readthedocs.io/en/release0.5.5/modes.html#fixed-rate-mode).
"""
mutable struct CompressedArraySeq{T,Nx} <: AbstractCompArraySeq
    data::Vector{UInt8}
    headpositions::Vector{Int64}
    tailpositions::Vector{Int64}
    spacedim::NTuple{Nx,Int32}
    timedim::Int32
    eltype::Type{T}
    tol::Float32
    precision::Int64
    rate::Int64

    function CompressedArraySeq(dtype::DataType, spacedim::Integer...; nt::Integer=1, rate::Int=0, tol::Real=0, precision::Int=0)
        data = Vector{UInt8}()
        headpositions = Int64[0] # trick to avoid checking for the first iteration in the append! function
        tailpositions = Int64[0] # which means the timedim = length(tailpositions) - 1
        eltype = dtype
        timedim = 0

        compSeq = new{dtype, length(spacedim)}(data, headpositions, tailpositions, spacedim, timedim, eltype, tol, precision, rate)

        expectedSize = length(zfp_compress(rand(dtype, spacedim...), write_header=false, rate=rate, tol=tol, precision=precision))
        if rate > 0 || precision > 0
            sizehint!(compSeq.data, nt*expectedSize+1)
        else tol > 0
            sizehint!(compSeq.data, Int64(round(nt*expectedSize*0.7)))
        end

        return compSeq
    end
end

"""
    size(compArray::AbstractCompArraySeq)

Returns the dimensions of the uncompressed array, with the last dimension being the time dimension.
"""
Base.size(compArray::AbstractCompArraySeq) = (compArray.spacedim..., compArray.timedim)

"""
    ndims(compArray::AbstractCompArraySeq)

Returns the number of dimensions of the uncompressed array, including the time dimension.
"""
Base.ndims(compArray::AbstractCompArraySeq) = length(compArray.spacedim) + 1

Base.IndexStyle(::Type{<:AbstractCompArraySeq}) = IndexLinear()

"""
    getindex(compArray::AbstractCompArraySeq, timeidx::Int)

Retrieve and decompress a single time slice from `compArray` at `timeidx`.
"""
Base.@propagate_inbounds function Base.getindex(compArray::CompressedArraySeq, timeidx::Int)
    @boundscheck timeidx <= compArray.timedim || throw(BoundsError(compArray, timeidx))

    decompArray = zeros(compArray.eltype, compArray.spacedim...)
    @inbounds zfp_decompress!(decompArray,
                              @views compArray.data[compArray.tailpositions[timeidx+1]:compArray.headpositions[timeidx+1]];
                              tol=compArray.tol, precision=compArray.precision, rate=compArray.rate)
    return decompArray
end

"""
    append!(compArray::CompressedArraySeq{T,N}, array::AbstractArray{T,N})

Append a new time slice to compArray, compressing array in the process.

# Arguments

    compArray::CompressedArraySeq{T,N}: Existing compressed array.
    array::AbstractArray{T,N}: Uncompressed array to append.
```
"""
function Base.append!(compArray::CompressedArraySeq{T,N}, array::AbstractArray{T,N}) where {T<:AbstractFloat, N}
    data = zfp_compress(array, write_header=false,
                        tol=compArray.tol, precision=compArray.precision, rate=compArray.rate,
                        nthreads=Threads.nthreads())
    fileSize = length(data)
    append!(compArray.data, data)
    push!(compArray.tailpositions, compArray.headpositions[end]+1)
    push!(compArray.headpositions, compArray.headpositions[end]+fileSize)
    compArray.timedim += 1
    return nothing
end

# function Base.read!(compArray::CompressedArraySeq{T,N}, outArray::AbstractArray{T,N}, timeidx::Int) where {T<:AbstractFloat, N}
    # @boundscheck timeidx <= compArray.timedim || throw(BoundsError(compArray, timeidx))

    # @inbounds zfp_decompress!(decompArray,
                              # compArray.data[compArray.tailpositions[timeidx+1]:compArray.headpositions[timeidx+1]];
                              # tol=compArray.tol, precision=compArray.precision, rate=compArray.rate)

    # return nothing
# end
