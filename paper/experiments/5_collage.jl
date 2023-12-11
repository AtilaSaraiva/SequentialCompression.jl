using CairoMakie
using SequentialCompression
using JSON

function add_axis_inset(pos=fig[1, 1]; backgroundcolor=:snow2,
    halign, valign, width=Relative(0.5),height=Relative(0.35),
    alignmode=Mixed(left=2, right=10), limits, axisOptions...)

    inset_box = Axis(pos; width, height, halign, valign, alignmode,
        backgroundcolor=backgroundcolor, axisOptions..., limits=limits)
    # bring content upfront
    translate!(inset_box.scene, 0, 0, 10)
    return inset_box
end

function compressData(Aoriginal,
                 tol1::Real=0,
                 tol2::Real=0,
                 inmemory::Bool=false)
    Acomp1 = SeqCompressor(Float64, size(Aoriginal)[1:2]..., inmemory=inmemory, tol=tol1)
    Acomp2 = SeqCompressor(Float64, size(Aoriginal)[1:2]..., inmemory=inmemory, tol=tol2)

    ny, nx, nt = size(Aoriginal)

    for it = 1:nt
        append!(Acomp1, Aoriginal[it])
        append!(Acomp2, Aoriginal[it])
    end
    return Acomp1, Acomp2
end

function collage(Aoriginal::SequentialCompression.AbstractCompArraySeq,
                 Acomp1::SequentialCompression.AbstractCompArraySeq,
                 Acomp2::SequentialCompression.AbstractCompArraySeq,
                 geometry::Dict,
                 numberOfSnapshots::Int=3;
                 tol1::Real=0,
                 tol2::Real=0,
                 boxes::Tuple=())

    @assert length(boxes) == numberOfSnapshots

    ny, nx, nt = size(Aoriginal)

    N = numberOfSnapshots
    M = 3

    cm_to_pt(cm) = cm .* 28.3465
    size_in_cm = (20, 18)
    size_in_pt = cm_to_pt(size_in_cm)
    fig = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), size=size_in_pt, fontsize=10)

    it = range(1//4*nt, 1//2*nt|>round|>Int, length=N)
    axisOptions = (:yreversed => true,
                   :xaxisposition => :top,
                   :xminorticksvisible => true,
                   :yminorticksvisible => true,
                   :ytickformat => x -> string.(round.(geometry["dy"] .* x, digits=2)),
                   :xtickformat => x -> string.(round.(geometry["dx"] .* x, digits=2)),
                   :titlesize => 14,
                  )

    traceAxis = Vector{Axis}(undef, N)

    for (n, it) in zip(1:N, it)
        it = round(it)
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
            Label(gl[1, 1], "time=$(time) s", padding = (5, 5, 5, 5))
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

        hideydecorations!(ax2, ticks=false, minorticks=false, grid=true)
        hideydecorations!(ax3, ticks=false, minorticks=false, grid=true)

        img1 = image!(ax1, Aoriginal[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)
        image!(ax2, Acomp1[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)
        image!(ax3, Acomp2[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)

        # Insets
        if boxes[n] != false
            ax1Inset = add_axis_inset(fig[n,1]; boxes[n]..., axisOptions...)
            ax2Inset = add_axis_inset(fig[n,2]; boxes[n]..., axisOptions...)
            ax3Inset = add_axis_inset(fig[n,3]; boxes[n]..., axisOptions...)
            image!(ax1Inset, Aoriginal[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)
            image!(ax2Inset, Acomp1[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)
            image!(ax3Inset, Acomp2[Int(it)]; colorrange=(vmin,vmax), colormap=:seismic)
        end

        Colorbar(fig[n,4], img1)

        colors = [:cyan, :magenta, :yellow]
        lines!(ax1, [nx÷2, nx÷2], [1, ny], color=colors[1], linewidth=1.5, linestyle=:dash)
        lines!(ax2, [nx÷2, nx÷2], [1, ny], color=colors[2], linewidth=1.5, linestyle=:dash)
        lines!(ax3, [nx÷2, nx÷2], [1, ny], color=colors[3], linewidth=1.5, linestyle=:dash)

        ax4 = Axis(fig[n,5], yreversed = true,
                   xaxisposition = :top,
                   ytickformat = x -> string.(round.(geometry["dy"] .* x, digits=2)),
                   yautolimitmargin=(0,0), xautolimitmargin=(0,0))

        traceAxis[n] = ax4

        if n == 1
            ax4.xlabel = "Amplitude"
        end

        linkyaxes!(ax1, ax4)

        lines!(ax4, Aoriginal[Int(it)][:, nx÷2], 1:ny, color=colors[1], linewidth=2, linestyle=:solid, label="Original")
        lines!(ax4, Acomp1[Int(it)][:, nx÷2], 1:ny, color=colors[2], linewidth=2, linestyle=:solid, label="tol = $(tol1)")
        lines!(ax4, Acomp2[Int(it)][:, nx÷2], 1:ny, color=colors[3], linewidth=2, linestyle=:dash, label="tol = $(tol2)")

    end

    colgap!(fig.layout, 2)
    rowgap!(fig.layout, 2)
    colgap!(fig.layout, 3,  5)
    colgap!(fig.layout, 4,  7)
    colsize!(fig.layout, 1, Relative(0.27))
    colsize!(fig.layout, 2, Relative(0.27))
    colsize!(fig.layout, 3, Relative(0.27))

    linkxaxes!(traceAxis...)

    CairoMakie.save("../figs/collage.pdf", fig, pt_per_unit=1)

    return nothing
end

Aoriginal = load("./wavefield.szfp")

geometry = open("./geometry.json") do file
    JSON.parse(file)
end

idx_x(x) = x / geometry["dx"]
idx_y(y) = y / geometry["dy"]

boxes = (
            false,
            Dict(
                :width=>Relative(0.9),
                :height=>Relative(0.35),
                :halign=>0.4,
                :valign=>0.3,
                :limits=>(idx_x(20), idx_x(30), idx_y(3), idx_y(9)),
                :xticks=>LinearTicks(2),
                :yticks=>LinearTicks(2),
                :yticklabelcolor=>:yellow,
                :xticklabelcolor=>:yellow,
            ),
            false
        )


tols = (1e-4, 1e-2)
Acomp1, Acomp2 = compressData(Aoriginal, tols...)
collage(Aoriginal, Acomp1, Acomp2, geometry, tol1=tols[1], tol2=tols[2], boxes=boxes)
