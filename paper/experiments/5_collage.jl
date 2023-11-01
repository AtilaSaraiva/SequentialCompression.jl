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

    ny, nx, nt = size(Aoriginal)

    for it = 1:nt
        append!(Acomp1, Aoriginal[it])
        append!(Acomp2, Aoriginal[it])
    end

    N = numberOfSnapshots
    M = 3

    size_in_inches = (5, 5)
    dpi = 600
    size_in_pixels = size_in_inches .* dpi
    fig = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), resolution=size_in_pixels)

    it = range(1, nt, length=N)
    axisOptions = (:yreversed => true,
                   :xaxisposition => :top,
                   :xminorticksvisible => true,
                   :yminorticksvisible => true,
                   :ytickformat => x -> string.(round.(geometry["dy"] .* x, digits=2)),
                   :xtickformat => x -> string.(round.(geometry["dx"] .* x, digits=2)),
                  )

    for (n, it) in zip(1:N, it)
        vmax = maximum(Aoriginal[Int(it)])
        vmin = minimum(Aoriginal[Int(it)])

        time = geometry["dt"] * (it - 1)
        time = round(time, digits=2)

        ax1 = Axis(fig[n,1]; axisOptions..., ylabel="Depth (m)")
        ax2 = Axis(fig[n,2]; axisOptions...)
        ax3 = Axis(fig[n,3]; axisOptions...)

        # This is a hack to get the time label in the bottom right corner
        for m = 1:3
            gl = GridLayout(fig[n, m], tellwidth = false, tellheight = false, halign = :right, valign = :bottom)
            Box(gl[1, 1], color = :bisque, strokecolor = :black, strokewidth = 2)
            Label(gl[1, 1], "time=$(time)", padding = (10, 10, 10, 10))
        end

        if n == 1
            ax1.title = "Original"
            ax2.title = "tol = $(tol1)"
            ax3.title = "tol = $(tol2)"
            for ax in (ax1, ax2, ax3)
                ax.xlabel = "Offset (m)"
            end
        else
            for ax in (ax1, ax2, ax3)
                hidexdecorations!(ax, ticks=false, minorticks=false, grid=false)
            end
        end

        hideydecorations!(ax2, ticks=false, grid=false)
        hideydecorations!(ax3, ticks=false, grid=false)

        img1 = image!(ax1, Aoriginal[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)
        image!(ax2, Acomp1[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)
        image!(ax3, Acomp2[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)

        Colorbar(fig[n,4], img1)

        colors = [:cyan, :magenta, :yellow]
        lines!(ax1, [nx÷2, nx÷2], [1, ny], color=colors[1], linewidth=1.5, linestyle=:dash)
        lines!(ax2, [nx÷2, nx÷2], [1, ny], color=colors[2], linewidth=1.5, linestyle=:dash)
        lines!(ax3, [nx÷2, nx÷2], [1, ny], color=colors[3], linewidth=1.5, linestyle=:dash)

        ax4 = Axis(fig[n,5], ylabel="Depth (m)", xlabel="Amplitude", yreversed = true,
                   xaxisposition = :top,
                   ytickformat = x -> string.(round.(geometry["dy"] .* x, digits=2)),
                   yautolimitmargin=(0,0), xautolimitmargin=(0,0))

        lines!(ax4, Aoriginal[Int(it)][:, nx÷2], 1:ny, color=colors[1], linewidth=3, linestyle=:solid, label="Original")
        lines!(ax4, Acomp1[Int(it)][:, nx÷2], 1:ny, color=colors[2], linewidth=3, linestyle=:solid, label="tol = $(tol1)")
        lines!(ax4, Acomp2[Int(it)][:, nx÷2], 1:ny, color=colors[3], linewidth=3, linestyle=:dash, label="tol = $(tol2)")
    end

    display(fig)

    return nothing
end

Aoriginal = load("./wavefield.szfp")

geometry = open("./geometry.json") do file
    JSON.parse(file)
end

collage(Aoriginal, geometry, tol1=1e-4, tol2=1e-2)
