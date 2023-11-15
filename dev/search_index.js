var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SequentialCompression","category":"page"},{"location":"#SequentialCompression","page":"Home","title":"SequentialCompression","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SequentialCompression.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SequentialCompression]","category":"page"},{"location":"#SequentialCompression.CompressedArraySeq","page":"Home","title":"SequentialCompression.CompressedArraySeq","text":"CompressedArraySeq{T,Nx}\n\nA mutable structure for storing time-dependent arrays in a compressed format.\n\nFields\n\ndata::Vector{UInt8}: Compressed data in byte form.\nheadpositions::Vector{Int64}: Positions of the beginning of each time slice in data.\ntailpositions::Vector{Int64}: Positions of the end of each time slice in data.\nspacedim::NTuple{Nx,Int32}: Dimensions of the spatial grid.\ntimedim::Int32: Number of time steps.\neltype::Type{T}: Element type of the uncompressed array.\ntol::Float32: Mean absolute error that is tolerated.\nprecision::Float32: Controls the precision, bounding a weak relative error.\nrate::Int64: Fixes the bits used per value.\n\n\n\n\n\n","category":"type"},{"location":"#SequentialCompression.CompressedMultiFileArraySeq","page":"Home","title":"SequentialCompression.CompressedMultiFileArraySeq","text":"CompressedMultiFileArraySeq{T,Nx}\n\nA compressed time-dependent array that is stored in multiple files, one per thread.\n\nFields\n\nfiles::Vector{IOStream}: IO object for each array slice.\nheadpositions::Vector{Int64}: Positions of the beginning of each time slice in data.\ntailpositions::Vector{Int64}: Positions of the end of each time slice in data.\nspacedim::NTuple{Nx,Int32}: Dimensions of the spatial grid.\ntimedim::Int32: Number of time steps.\neltype::Type{T}: Element type of the uncompressed array.\ntol::Float32: Mean absolute error that is tolerated.\nprecision::Float32: Controls the precision, bounding a weak relative error.\nrate::Int64: Fixes the bits used per value.\n\nArguments exclusive for the constructor\n\nfilepaths::Union{Vector{String}, String}=\"/tmp/seqcomp\": Path(s) to the files where the compressed data will be stored. If only one string is passed, the same path will be used for all threads.\n\n\n\n\n\n","category":"type"},{"location":"#Base.append!-Union{Tuple{N}, Tuple{T}, Tuple{SequentialCompression.CompressedArraySeq{T, N}, AbstractArray{T, N}}} where {T<:AbstractFloat, N}","page":"Home","title":"Base.append!","text":"append!(compArray::CompressedArraySeq{T,N}, array::AbstractArray{T,N})\n\nAppend a new time slice to compArray, compressing array in the process.\n\nArguments\n\ncompArray::CompressedArraySeq{T,N}: Existing compressed array.\narray::AbstractArray{T,N}: Uncompressed array to append.\n\n```\n\n\n\n\n\n","category":"method"},{"location":"#Base.getindex-Tuple{SequentialCompression.CompressedArraySeq, Int64}","page":"Home","title":"Base.getindex","text":"getindex(compArray::AbstractCompArraySeq, timeidx::Int)\n\nRetrieve and decompress a single time slice from compArray at timeidx.\n\n\n\n\n\n","category":"method"},{"location":"#Base.ndims-Tuple{SequentialCompression.AbstractCompArraySeq}","page":"Home","title":"Base.ndims","text":"ndims(compArray::AbstractCompArraySeq)\n\nReturns the number of dimensions of the uncompressed array, including the time dimension.\n\n\n\n\n\n","category":"method"},{"location":"#Base.size-Tuple{SequentialCompression.AbstractCompArraySeq}","page":"Home","title":"Base.size","text":"size(compArray::AbstractCompArraySeq)\n\nReturns the dimensions of the uncompressed array, with the last dimension being the time dimension.\n\n\n\n\n\n","category":"method"},{"location":"#SequentialCompression.SeqCompressor-Tuple{DataType, Vararg{Integer}}","page":"Home","title":"SequentialCompression.SeqCompressor","text":"SeqCompressor(dtype::DataType, spacedim::Integer...;\n              rate::Int=0, tol::Real=0, precision::Real=0,\n              filepaths::Union{Vector{String}, String}=\"\",\n              envVarPath::String=\"\")\n\nConstruct a CompressedArraySeq or CompressedMultiFileArraySeq depending on the arguments.\n\nArguments\n\ndtype::DataType: the type of the array to be compressed\nspacedim::Integer...: the dimensions of the array to be compressed\ninmemory::Bool=true: whether the compressed data will be stored in memory or in disk\nrate::Int64: Fixes the bits used per value.\ntol::Float32: Mean absolute error that is tolerated.\nprecision::Float32: Controls the precision, bounding a weak relative error.\nfilepaths::Union{Vector{String}, String}=\"\": the path(s) to the files to be compressed\nenvVarPath::String=\"\": the name of the environment variable that contains the path to the files to be compressed\n\nYou have the option of passing an environment variable, a file path, a vector of file paths, or nothing. If you pass a vector of file paths, the number of paths must be equal to the number of threads. If you pass a single file path, the same path will be used for all threads. If you pass an environment variable, the file path will be extracted from it. It might be useful if you are using a SLURM job scheduler, for example, since the local disk of the node can be accessed by ENV[\"SLURM_TMPDIR\"].\n\nExample\n\njulia> using SequentialCompression\n\njulia> A = SeqCompressor(Float64, 4, 4)\nSequentialCompression.CompressedArraySeq{Float64, 2}(UInt8[], [0], [0], (4, 4), 0, Float64, 0.0f0, 0, 0)\n\njulia> A.timedim\n0\n\njulia> size(A)\n(4, 4, 0)\n\njulia> append!(A, ones(Float64, 4, 4));\n\njulia> A[1]\n4×4 Matrix{Float64}:\n 1.0  1.0  1.0  1.0\n 1.0  1.0  1.0  1.0\n 1.0  1.0  1.0  1.0\n 1.0  1.0  1.0  1.0\n\njulia> size(A)\n(4, 4, 1)\n\n\n\n\n\n","category":"method"}]
}
