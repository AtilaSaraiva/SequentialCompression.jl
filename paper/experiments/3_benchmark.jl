using SequentialCompression
using ZfpCompression
using CairoMakie
using DataFrames

const N1, N2 = 2000, 2000

function compressionthroughput(; rate::Int=0, tol::Real=0, precision::Int=0, inmemory::Bool=true)
    dummyData = rand(N1, N2,10)
    if inmemory
        cp = SeqCompressor(Float64, N1, N2, rate=rate, tol=tol, precision=precision, inmemory=inmemory)
    else
        cp = SeqCompressor(Float64, N1, N2, rate=rate, tol=tol, precision=precision, inmemory=inmemory, filepaths="./seqcomp")
    end

    wtime0 = Base.time()
    for it=1:size(dummyData, 3)
        append!(cp, dummyData[:,:,it])
    end
    wtime    = Base.time()-wtime0

    throughput = sizeof(dummyData) * 1e-6 / wtime
    return throughput, cp
end

function decompressionthroughput(cp)
    wtime0 = Base.time()
    for it=1:size(cp)[3]
        cp[it]
    end
    wtime    = Base.time()-wtime0

    throughput = *(size(cp)...)*sizeof(cp.eltype) * 1e-6 / wtime
    return throughput
end

function compressionthroughput2()
    dummyData = rand(1000,1000,10)
    cp = Vector{Vector{UInt8}}(undef, 10)
    wtime0 = Base.time()
    for it=1:size(dummyData, 3)
        cp[it] = zfp_compress(dummyData[:,:,it], write_header=false, nthreads=8)
    end
    wtime    = Base.time()-wtime0

    throughput = sizeof(dummyData) * 1e-6 / wtime
    return throughput, cp
end

function decompressionthroughput2(cp)
    wtime0 = Base.time()
    decompData = zeros(1000,1000,10)
    for it=1:length(cp)
        A = zeros(1000,1000)
        zfp_decompress!(A, cp[it])
        decompData[:,:,it] .= A
    end
    wtime    = Base.time()-wtime0

    throughput = sizeof(decompData) * 1e-6 / wtime
    return throughput
end

function compressionthroughput3(; rate::Int=0, tol::Real=0, precision::Int=0, inmemory::Bool=true)
    dummyData = rand(1000,1000,10)
    cp = Vector{Vector{UInt8}}(undef, 10)
    wtime0 = Base.time()
    for it=1:size(dummyData, 3)
        cp[it] = zfp_compress(dummyData[:,:,it], rate=rate, tol=tol, precision=precision)
    end
    wtime    = Base.time()-wtime0

    throughput = sizeof(dummyData) * 1e-6 / wtime
    return throughput, cp
end

function decompressionthroughput3(cp)
    decompData = zeros(1000,1000)
    wtime0 = Base.time()
    for it=1:10
        decompData[:,:] .= zfp_decompress(cp[it])
    end
    wtime    = Base.time()-wtime0

    throughput = sizeof(decompData) * 1e-6 / wtime
    return throughput
end

function throughputs(; rate::Int=0, tol::Real=0, precision::Int=0, inmemory::Bool=true)
    comp_throughput, cp = compressionthroughput(rate=rate, tol=tol, precision=precision, inmemory=inmemory)
    decomp_throughput = decompressionthroughput(cp)
    return comp_throughput, decomp_throughput
end

function rateTest()
    rates = 4:8:52
    df = DataFrame(rate=Int64[], media=String[], compThroughput=Float64[], decompThroughput=Float64[])

    for (i,rate) in enumerate(rates)
        comp_throughput, decomp_throughput = throughputs(rate=rate, inmemory=true)
        push!(df, (rate, "inmemory", comp_throughput, decomp_throughput))
    end

    for (i,rate) in enumerate(rates)
        comp_throughput, decomp_throughput = throughputs(rate=rate, inmemory=false)
        push!(df, (rate, "disk", comp_throughput, decomp_throughput))
    end

    gdf = groupby(df, Cols(:rate, :media))
    return combine(gdf, :compThroughput => maximum, :decompThroughput => maximum, renamecols=false)
end

function precisionTest()
    precisions = 4:8:52
    df = DataFrame(precision=Int64[], media=String[], compThroughput=Float64[], decompThroughput=Float64[])

    for (i,precision) in enumerate(precisions)
        comp_throughput, decomp_throughput = throughputs(precision=precision, inmemory=true)
        push!(df, (precision, "inmemory", comp_throughput, decomp_throughput))
    end

    for (i,precision) in enumerate(precisions)
        comp_throughput, decomp_throughput = throughputs(precision=precision, inmemory=false)
        push!(df, (precision, "disk", comp_throughput, decomp_throughput))
    end

    gdf = groupby(df, Cols(:precision, :media))
    return combine(gdf, :compThroughput => maximum, :decompThroughput => maximum, renamecols=false)
end

function tolTest()
    tols = [ 10.0^(-i) for i=1:1:10 ]
    df = DataFrame(tol=Float64[], media=String[], compThroughput=Float64[], decompThroughput=Float64[])

    for (i,tol) in enumerate(tols)
        comp_throughput, decomp_throughput = throughputs(tol=tol, inmemory=true)
        push!(df, (tol, "inmemory", comp_throughput, decomp_throughput))
    end

    for (i,tol) in enumerate(tols)
        comp_throughput, decomp_throughput = throughputs(tol=tol, inmemory=false)
        push!(df, (tol, "disk", comp_throughput, decomp_throughput))
    end

    gdf = groupby(df, Cols(:tol, :media))
    return combine(gdf, :compThroughput => maximum, :decompThroughput => maximum, renamecols=false)
end

function plot(df_rate, df_tol, df_precision)
    cm_to_pt(cm) = cm .* 28.3465
    size_in_cm = (18, 15)
    size_in_pt = cm_to_pt(size_in_cm)
    fig = Figure(size=size_in_pt)
    ax1 = Axis(fig[1,1],
               xlabel="rate (bps)",
               ylabel="throughput (MiB/s)",
               title="throughput x rate",
               xautolimitmargin=(0,0),
               xticks = (4:4:52)
             )
    ax2 = Axis(fig[2,1],
               xlabel="tolerance",
               ylabel="throughput (MiB/s)",
               title="throughput x tolerance",
               xscale=log10,
               xautolimitmargin=(0,0),
               xminorticksvisible=false,
               xreversed=true
             )
    ax3 = Axis(fig[3,1],
               xlabel="precision (bit planes)",
               ylabel="throughput (MiB/s)",
               title="throughput x precision",
               xautolimitmargin=(0,0),
               xticks = (4:4:52)
             )


    colors = [:blue, :red, :purple, :orange]
    for (df, ax, var) in zip((df_rate, df_tol, df_precision), (ax1, ax2, ax3), (:rate, :tol, :precision))
        let df_inmemory = filter(:media => x -> x=="inmemory", df),
            df_disk = filter(:media => x -> x=="disk", df),
            ln1 = lines!(ax, df_inmemory[!,var], df_inmemory.compThroughput, color=colors[1])
            ln2 = lines!(ax, df_inmemory[!,var], df_inmemory.decompThroughput, color=colors[2])
            ln3 = lines!(ax, df_disk[!,var], df_disk.compThroughput, color=colors[3])
            ln4 = lines!(ax, df_disk[!,var], df_disk.decompThroughput, color=colors[4])
        end
    end

    group1 = [LineElement(color=color, linestype=nothing) for color in colors[1:2]]
    group2 = [LineElement(color=color, linestype=nothing) for color in colors[3:4]]
    labels = ["Compression throughput", "Decompression throughput"]
    legend = Legend(fig, [group1, group2], [labels, labels], ["In-memory", "Disk"], framevisible=false)
    fig[1:3,2] = legend

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 10)
    return fig
end

publication_theme() = Theme(
    fontsize=10,
    font="CMU Serif",
    figure_padding=8,
    Axis=(
        xgridstyle=:dash, ygridstyle=:dash,
        xminorticksvisible=true,
        ),
    Legend=(framecolor=(:black, 0.5), backgroundcolor=(:white, 0.5), framevisible=false, tellheight=true, tellwidth=false, labelsize=8)
)

precision_df = precisionTest()
tol_df = tolTest()
rate_df = rateTest()
fig = with_theme(publication_theme()) do
    plot(rate_df, tol_df, precision_df)
end

CairoMakie.save("../figs/throughput.pdf", fig, pt_per_unit=1)
