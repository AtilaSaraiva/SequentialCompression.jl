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


function plot()
    cm_to_pt(cm) = cm .* 28.3465
    size_in_cm = (8.16, 13)
    size_in_pt = cm_to_pt(size_in_cm)
    fig = Figure(size=size_in_pt)
    ax1 = Axis(fig[1,1],
              xlabel = "Finite difference iterations",
              ylabel = "Size of compressed array (Mb)"
             )

    ln1 = lines!(ax1, sizes_lossless.wave)
    ln2 = lines!(ax1, sizes_lossless.random)
    ln3 = lines!(ax1, sizes_lossy.wave)
    ln4 = lines!(ax1, sizes_lossy.random)

    group1 = [LineElement(color=ln.color, linestype=ln.linestyle) for ln in (ln1, ln2)]
    group2 = [LineElement(color=ln.color, linestype=ln.linestyle) for ln in (ln3, ln4)]
    labels = ["wavefield compressed snapshot size", "random arrays compressed snapshot size"]
    legend = Legend(fig[2,1], [group1, group2], [labels, labels], ["Lossless", "Lossy"], framevisible=false)

    return fig
end

publication_theme = Theme(
    fontsize=10,
    font="CMU Serif",
    figure_padding=8,
    Axis=(
        xgridstyle=:dash, ygridstyle=:dash,
        xminorticksvisible=true,
        ),
    Legend=(framecolor=(:black, 0.5), backgroundcolor=(:white, 0.5), framevisible=false, tellheight=true, tellwidth=false, labelsize=8, orientation=:vertical, titlesize=9),
)

fig = with_theme(plot, publication_theme)
CairoMakie.save("../figs/compressedSizeAnalysis.pdf", fig)
