using SequentialCompression
using CairoMakie
import Base

function compressData(Aoriginal, rate::Int=1)
    Ainmemory  = SeqCompressor(Float64, size(Aoriginal)[1:2]..., inmemory = true, rate = rate)
    Amultifile = SeqCompressor(Float64, size(Aoriginal)[1:2]..., inmemory = false, rate = rate)

    ny, nx, nt = size(Aoriginal)

    for it = 1:nt
        append!(Ainmemory, Aoriginal[it])
        append!(Amultifile, Aoriginal[it])
    end

    return Ainmemory, Amultifile
end


Aoriginal = load("./wavefield.szfp")
shape = size(Aoriginal)

rates = 1:9:64

multifileSizes = zeros(Int32, length(rates))
inmemorySizes = zeros(Int32, length(rates))
for (i, rate) in enumerate(rates)
    Ainmemory, Amultifile = compressData(Aoriginal, rate)
    map(flush, Amultifile.files)
    filesTotalSize = mapreduce(filesize, +, Amultifile.files)
    inmemoryStructSize, multifileStructSize = map(Base.summarysize, (Ainmemory, Amultifile))
    multifileStructSize = filesTotalSize + multifileStructSize
    multifileSizes[i] = multifileStructSize
    inmemorySizes[i] = inmemoryStructSize
end

originalFileSize = *(size(Aoriginal)...) * sizeof(Aoriginal.eltype)
compressionRatio(compressedSize) = originalFileSize / compressedSize
spaceSaving(compressedSize) = 1 - compressedSize / originalFileSize

publication_theme() = Theme(
    fontsize=10,
    font="CMU Serif",
    figure_padding=8,
    Axis=(
        xgridstyle=:dash, ygridstyle=:dash,
        xtickalign=1, ytickalign=1,
        xminortickalign=1, yminortickalign=1,
        xlabel="rate (bpd)", ylabel="space saving",
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(5),
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(5),
        xminortickcolor=:gray,
        yminortickcolor=:gray,
        ),
    ScatterLines = (markersize=6,),
    Legend=(framecolor=(:black, 0.5), backgroundcolor=(:white, 0.5)),
    Colorbar=(ticksize=16, tickalign=1, spinewidth=0.5),
)

function plot()
    cm_to_pt(cm) = cm .* 28.3465
    size_in_cm = (8.16, 8.16)
    size_in_pt = cm_to_pt(size_in_cm)
    fig = Figure(size=size_in_pt)
    ax = Axis(fig[1,1])
    scatterlines!(ax, rates, spaceSaving.(multifileSizes), label="Multifile")
    scatterlines!(ax, rates, spaceSaving.(inmemorySizes), label="Inmemory")
    axislegend(;position=:lb, rowgap=1)
    return fig
end


fig = with_theme(plot, publication_theme())
CairoMakie.save("../figs/ratio.pdf", fig, pt_per_unit=1)
