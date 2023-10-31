using CairoMakie
using SequentialCompression
using JSON


function collage(Aoriginal::SequentialCompression.AbstractCompArraySeq,
                 geometry::Dict,
                 numberOfSnapshots::Int=4;
                 tol1::Real=0,
                 tol2::Real=0,
                 inmemory::Bool=true)

    Acomp1 = SeqCompressor(Float64, size(Aoriginal)[1:2]..., inmemory=false, tol=tol1)
    Acomp2 = SeqCompressor(Float64, size(Aoriginal)[1:2]..., inmemory=false, tol=tol2)

    nt = size(Aoriginal)[end]

    for it = 1:nt
        append!(Acomp1, Aoriginal[it])
        append!(Acomp2, Aoriginal[it])
    end

    N = numberOfSnapshots
    M = 3

    size_in_inches = (4, 5)
    dpi = 300
    size_in_pixels = size_in_inches .* dpi
    fig = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), resolution=size_in_pixels)

    it = range(1, nt, length=N)
    axisOptions = (:yreversed => true,
                   :xaxisposition => :top)

    for (n, it) in zip(1:N, it)
        vmax = maximum(Aoriginal[Int(it)])
        vmin = minimum(Aoriginal[Int(it)])
        ax1 = Axis(fig[n,1]; axisOptions..., ylabel="Depth (m)")
        ax2 = Axis(fig[n,2]; axisOptions...)
        ax3 = Axis(fig[n,3]; axisOptions...)

        if n == 1
            ax1.title = "Original"
            ax2.title = "tol = $(tol1)"
            ax3.title = "tol = $(tol2)"
            for ax in (ax1, ax2, ax3)
                ax.xlabel = "Offset (m)"
            end
        end

        img1 = image!(ax1, Aoriginal[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)
        image!(ax2, Acomp1[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)
        image!(ax3, Acomp2[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)
        Colorbar(fig[n,4], img1)
    end

    display(fig)

    return nothing
end

Aoriginal = load("./wavefield.szfp")

geometry = open("./geometry.json") do file
    JSON.parse(file)
end

collage(Aoriginal, geometry, tol1=1e-2, tol2=1e-4)
