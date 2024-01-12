using DataFrames
using CSV
using CairoMakie
using Statistics
import Makie: pseudolog10, Symlog10

"""
This function filters the dataframe to grab only rows where columnToFilter is bigger than zero,
and then groups it by the colsToGroup, and then combine the groups while calculating the minimum and
maximum value of the wtime, wtime_it and T_eff
"""
function processData(df::DataFrame, colsToGroup::Cols, columnToFilter::Union{Symbol,Nothing}=nothing)
    if !isnothing(columnToFilter)
        filteredData = filter(columnToFilter => x -> x > 0, df)
    else
        filteredData = df
    end
    gdf = groupby(filteredData, colsToGroup)
    columnsToModify = [:wtime, :wtime_it, :T_eff]
    minMaxData = combine(gdf, columnsToModify .=> minimum, columnsToModify .=> maximum)
    return minMaxData
end


function plot()
    dfLimits = DataFrame(CSV.File("dfLimits.csv"))
    dfCompMultifile = DataFrame(CSV.File("dfCompMultifile.csv"))
    dfCompInmem = DataFrame(CSV.File("dfCompInmem.csv"))

    cm_to_pt(cm) = cm .* 28.3465
    size_in_cm = (8.16, 15)
    size_in_pt = cm_to_pt(size_in_cm)

    multifileData = processData(dfCompMultifile, Cols(:rate, :nthreads), :rate)
    inmemoryData = processData(dfCompInmem, Cols(:rate), :rate)
    comparissonData = processData(dfLimits, Cols(:media))

    maxThreads = maximum(multifileData.nthreads)
    multifileDataMaxThreads = filter(:nthreads => nthreads -> nthreads == maxThreads, multifileData)

    f = Figure(size=size_in_pt)

    # A
    ga = f[1,1] = GridLayout()

    ax1 = Axis(ga[1,1], ylabel="nthreads", xlabel="rate (bpd)")

    hm = heatmap!(ax1, multifileData.rate, multifileData.nthreads, 1000 .* multifileData.wtime_it_minimum)
    xlims!(ax1, minimum(multifileData.rate), maximum(multifileData.rate))
    ylims!(ax1, minimum(multifileData.nthreads), maximum(multifileData.nthreads))
    Colorbar(ga[1,2], hm; label = "exec. time per iter. (ms/it)", vertical=true)
    colgap!(ga, 8)

    # B
    gb = f[2,1] = GridLayout()
    ax2 = Axis(gb[1,1], alignmode=Mixed(right=0), ylabel="time / iter (ms/it)", xlabel="rate (bpd)")

    ln1 = lines!(ax2, multifileDataMaxThreads.rate, 1000 .* multifileDataMaxThreads.wtime_it_minimum, label="multifile comp.")
    band!(ax2, multifileDataMaxThreads.rate, 1000 .* multifileDataMaxThreads.wtime_it_minimum, 1000 .* multifileDataMaxThreads.wtime_it_maximum)

    ln2 = lines!(ax2, inmemoryData.rate, 1000 .* inmemoryData.wtime_it_minimum, label="in memory comp.")
    band!(ax2, inmemoryData.rate, 1000 .* inmemoryData.wtime_it_minimum, 1000 .* inmemoryData.wtime_it_maximum)

    let storage = filter(:media => media -> media == "disk", comparissonData),
        memory = filter(:media => media -> media == "memory", comparissonData)

        ln4 = lines!(ax2, inmemoryData.rate, fill(1000 * storage.wtime_it_minimum[1], length(inmemoryData.rate)), label="disk")

        ln5 = lines!(ax2, inmemoryData.rate, fill(1000 * memory.wtime_it_minimum[1], length(inmemoryData.rate)), label="array")

    end

    multifileData = processData(dfCompMultifile, Cols(:tol, :nthreads), :tol)
    inmemoryData = processData(dfCompInmem, Cols(:tol), :tol)

    maxThreads = maximum(multifileData.nthreads)
    multifileDataMaxThreads = filter(:nthreads => nthreads -> nthreads == maxThreads, multifileData)

    # # C

    # gc = f[1,2] = GridLayout()

    # ax3 = Axis(gc[1,1], ylabel="nthreads", xlabel="rate (bpd)", xscale=pseudolog10)
    # # xlims!(ax3, minimum(multifileData.tol), maximum(multifileData.tol))

    # hm = heatmap!(ax3, multifileData.tol, multifileData.nthreads, 1000 .* multifileData.wtime_it_minimum)
    # # Colorbar(gc[1,2], hm; label = "exec. time per iter. (ms/it)", vertical=true)
    # colgap!(gc, 8)

    # # D

    # gd = f[2,2] = GridLayout()
    # ax4 = Axis(gd[1,1], alignmode=Mixed(right=0), ylabel="time / iter (ms/it)", xlabel="tol (bpd)")

    # ln1 = lines!(ax4, multifileDataMaxThreads.tol, 1000 .* multifileDataMaxThreads.wtime_it_minimum, label="multifile comp.")
    # band!(ax2, multifileDataMaxThreads.tol, 1000 .* multifileDataMaxThreads.wtime_it_minimum, 1000 .* multifileDataMaxThreads.wtime_it_maximum)

    # ln2 = lines!(ax4, inmemoryData.tol, 1000 .* inmemoryData.wtime_it_minimum, label="in memory comp.")
    # band!(ax2, inmemoryData.tol, 1000 .* inmemoryData.wtime_it_minimum, 1000 .* inmemoryData.wtime_it_maximum)

    # let storage = filter(:media => media -> media == "disk", comparissonData),
        # memory = filter(:media => media -> media == "memory", comparissonData)

        # ln4 = lines!(ax4, inmemoryData.tol, fill(1000 * storage.wtime_it_minimum[1], length(inmemoryData.tol)), label="disk")

        # ln5 = lines!(ax4, inmemoryData.tol, fill(1000 * memory.wtime_it_minimum[1], length(inmemoryData.tol)), label="array")

    # end


    # rowsize!(f.layout, 1, Relative(0.5))

    rowgap!(f.layout, 10)

    Legend(f[3,1], ax2)

    label = [L"\textbf{(a)}", L"\textbf{(b)}", L"\textbf{(c)}", L"\textbf{(d)}"]
    Label(ga[1,1:2, Bottom()], label[1],
          fontsize = 12,
          padding = (0,0, 0,40),
          halign = :center)

    Label(gb[1,1, Bottom()], label[2],
          fontsize = 12,
          padding = (0,0, 0,40),
          halign = :center)

    return f
end

publication_theme() = Theme(
    fontsize=10,
    font="CMU Serif",
    figure_padding=8,
    Axis=(
        xgridstyle=:dash, ygridstyle=:dash,
        xminorticksvisible=true,
        ),
    Legend=(framecolor=(:black, 0.5), backgroundcolor=(:white, 0.5), framevisible=false, tellheight=true, tellwidth=false, labelsize=8,
            nbanks=2, orientation=:horizontal),
)


fig = with_theme(plot, publication_theme())

CairoMakie.save("../figs/speeddiff.pdf", fig, pt_per_unit=1)
