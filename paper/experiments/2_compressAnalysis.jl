using SequentialCompression
using ZfpCompression
using Random
using CairoMakie

struct Sizes
    wave::Vector{Float64}
    random::Vector{Float64}
end

function sizeComp(compress::Function)
    Ac = load("wavefield.szfp")

    sizes = Float64[]
    for i=1:size(Ac)[3]
       push!(sizes, compress(Ac[i]) |> sizeof)
    end

    Random.seed!(1234)

    dims = size(Ac)[1:2]

    create_random_arrays(dims) = rand(dims...)

    sizes_random = Float64[]
    for i=1:size(Ac)[3]
        push!(sizes_random, compress(create_random_arrays(dims)) |> sizeof)
    end

    sizes = map(x -> x/1024^2, sizes)
    sizes_random = map(x -> x/1024^2, sizes_random)

    return Sizes(sizes, sizes_random)
end

sizes_lossless = sizeComp(A -> zfp_compress(A))
sizes_lossy = sizeComp(A -> zfp_compress(A, tol=1e-5))

labels = ["sizes for each wavefield snapshot", "sizes for random arrays"]

fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98))
ax1 = Axis(fig[1:2,1],
          title = "Compressed array size: wavefield vs random arrays",
          xlabel = "Finite difference iterations",
          ylabel = "Size of compressed array (Mb)"
         )

ln1 = lines!(ax1, sizes_lossless.wave)
ln2 = lines!(ax1, sizes_lossless.random)
ln3 = lines!(ax1, sizes_lossy.wave)
ln4 = lines!(ax1, sizes_lossy.random)
legend = Legend(fig, [[ln1, ln2], [ln3, ln4]], [labels, labels], ["Lossless", "Lossy"], framevisible=false)

fig[1:2,2] = legend

display(fig)
CairoMakie.save("../figs/compressedSizeAnalysis.pdf", fig)
