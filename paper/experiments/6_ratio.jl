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

function compressData(Aoriginal, tol::Float64=1)
    Ainmemory  = SeqCompressor(Float64, size(Aoriginal)[1:2]..., inmemory = true, tol = tol)
    Amultifile = SeqCompressor(Float64, size(Aoriginal)[1:2]..., inmemory = false, tol = tol)

    ny, nx, nt = size(Aoriginal)

    for it = 1:nt
        append!(Ainmemory, Aoriginal[it])
        append!(Amultifile, Aoriginal[it])
    end

    return Ainmemory, Amultifile
end

function test(Aoriginal::SequentialCompression.AbstractCompArraySeq, values::AbstractVector)
    multifileSizes = zeros(Int64, length(values))
    inmemorySizes = zeros(Int64, length(values))
    for (i, value) in enumerate(values)
        Ainmemory, Amultifile = compressData(Aoriginal, value)
        map(flush, Amultifile.files)
        filesTotalSize = mapreduce(filesize, +, Amultifile.files)
        inmemoryStructSize, multifileStructSize = map(Base.summarysize, (Ainmemory, Amultifile))
        multifileStructSize = filesTotalSize + multifileStructSize
        multifileSizes[i] = multifileStructSize
        inmemorySizes[i] = inmemoryStructSize
    end

    return inmemorySizes, multifileSizes
end

Aoriginal = load("./wavefield.szfp")
shape = size(Aoriginal)


originalFileSize = *(size(Aoriginal)...) * sizeof(Aoriginal.eltype)
compressionRatio(compressedSize) = originalFileSize / compressedSize
spaceSaving(compressedSize) = 1 - compressedSize / originalFileSize

rates = 4:8:52
tols = [ 10.0^(-i) for i=1:1:10 ]
# inmemorySizesRate, multifileSizesRate = test(Aoriginal, rates)
# inmemorySizesTol, multifileSizesTol = test(Aoriginal, tols)

publication_theme() = Theme(
    fontsize=10,
    font="CMU Serif",
    figure_padding=8,
    Axis=(
        xgridstyle=:dash, ygridstyle=:dash,
        xtickalign=1, ytickalign=1,
        xminortickalign=1, yminortickalign=1,
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
    size_in_cm = (8.16, 17)
    size_in_pt = cm_to_pt(size_in_cm)
    fig = Figure(size=size_in_pt)
    ax1 = Axis(fig[1,1],
               xlabel="rate (bps)",
               ylabel="throughput (MiB/s)",
               title="space saving x rate",
               xautolimitmargin=(0,0),
               xticks = (4:4:52)
             )
    ax2 = Axis(fig[2,1],
               xlabel="tolerance",
               ylabel="space saving",
               title="space saving x tolerance",
               xscale=log10,
               xautolimitmargin=(0,0),
               xminorticksvisible=false,
               xreversed=true
             )

    scatterlines!(ax1, rates, spaceSaving.(multifileSizesRate), label="Multifile")
    scatterlines!(ax1, rates, spaceSaving.(inmemorySizesRate), label="Inmemory")

    scatterlines!(ax2, tols, spaceSaving.(multifileSizesTol), label="Multifile")
    scatterlines!(ax2, tols, spaceSaving.(inmemorySizesTol), label="Inmemory")

    axislegend(ax1; position=:lb, rowgap=1)
    axislegend(ax2; position=:lb, rowgap=1)

    label = [L"\textbf{(a)}", L"\textbf{(b)}", L"\textbf{(c)}", L"\textbf{(d)}"]
    Label(fig[1,1, Bottom()], label[1],
          fontsize = 12,
          padding = (0,0, 0,40),
          halign = :center)

    Label(fig[2,1, Bottom()], label[2],
          fontsize = 12,
          padding = (0,0, 0,40),
          halign = :center)
    return fig
end


fig = with_theme(plot, publication_theme())
CairoMakie.save("../figs/ratio.pdf", fig, pt_per_unit=1)
