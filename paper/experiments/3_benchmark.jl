using SequentialCompression
using ZfpCompression
using CairoMakie
using DataFrames

function compressionthroughput(; rate::Int=0, tol::Real=0, precision::Int=0, inmemory::Bool=true)
    dummyData = rand(1000,1000,10)
    cp = SeqCompressor(Float64, 1000, 1000, rate=rate, tol=tol, precision=precision, inmemory=inmemory)

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
    decompData = cp[:]
    wtime    = Base.time()-wtime0

    throughput = sizeof(decompData) * 1e-6 / wtime
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
    wtime0 = Base.time()
    decompData = zeros(1000,1000,10)
    for it=1:length(cp)
        decompData[:,:,it] .= zfp_decompress(cp[it])
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
    rates = 1:1:64
    df = DataFrame()
    df.rate = rates |> collect

    values = zeros(length(df.rate),2)
    for (i,rate) in enumerate(rates)
        values[i,:] .= throughputs(rate=rate, inmemory=true)
    end

    df.compression_inmemory = values[:, 1]
    df.decompression_inmemory = values[:, 2]

    for (i,rate) in enumerate(rates)
        values[i,:] .= throughputs(rate=rate, inmemory=false)
    end

    df.compression_disk = values[:, 1]
    df.decompression_disk = values[:, 2]

    return df
end

function precisionTest()
    precisions = 1:1:64
    df = DataFrame()
    df.precision = precisions |> collect

    values = zeros(length(df.precision),2)
    for (i,precision) in enumerate(precisions)
        values[i,:] .= throughputs(precision=precision, inmemory=true)
    end

    df.compression_inmemory = values[:, 1]
    df.decompression_inmemory = values[:, 2]

    for (i,precision) in enumerate(precisions)
        values[i,:] .= throughputs(precision=precision, inmemory=true)
    end

    df.compression_disk = values[:, 1]
    df.decompression_disk = values[:, 2]

    return df
end

function tolTest()
    tols = [ 10.0^(-i) for i=1:10 ]
    df = DataFrame()
    df.tol = tols

    values = zeros(length(df.tol),2)
    for (i,tol) in enumerate(tols)
        values[i,:] .= throughputs(tol=tol, inmemory=true)
    end

    df.compression_inmemory = values[:, 1]
    df.decompression_inmemory = values[:, 2]

    for (i,tol) in enumerate(tols)
        values[i,:] .= throughputs(tol=tol, inmemory=false)
    end

    df.compression_disk = values[:, 1]
    df.decompression_disk = values[:, 2]

    return df
end

function plot(df_rate, df_tol, df_precision)
    fig = Figure()
    ax1 = Axis(fig[1,1],
               xlabel="rate (bps)",
               ylabel="throughput (MiB/s)",
               title="Throughput x rate",
               xautolimitmargin=(0,0),
               xticks = (1:4:64)
             )
    ax2 = Axis(fig[2,1],
               xlabel="tolerance",
               ylabel="throughput (MiB/s)",
               title="Throughput x tolerance",
               xscale=log10,
               xautolimitmargin=(0,0),
             )
    ax3 = Axis(fig[3,1],
               xlabel="precision (bit planes)",
               ylabel="throughput (MiB/s)",
               title="Throughput x precision",
               xautolimitmargin=(0,0),
               xticks = (1:4:64)
             )

    colors = [:blue, :red, :green, :orange]
    ln1 = lines!(ax1, df_rate.rate, df_rate.compression_inmemory, color=colors[1])
    ln2 = lines!(ax1, df_rate.rate, df_rate.decompression_inmemory, color=colors[2])
    ln3 = lines!(ax1, df_rate.rate, df_rate.compression_disk, color=colors[3])
    ln4 = lines!(ax1, df_rate.rate, df_rate.decompression_disk, color=colors[4])

    lines!(ax2, df_tol.tol, df_tol.compression_inmemory, color=colors[1])
    lines!(ax2, df_tol.tol, df_tol.decompression_inmemory, color=colors[2])
    lines!(ax2, df_tol.tol, df_tol.compression_disk, color=colors[3])
    lines!(ax2, df_tol.tol, df_tol.decompression_disk, color=colors[4])

    lines!(ax3, df_precision.precision, df_precision.compression_inmemory, color=colors[1])
    lines!(ax3, df_precision.precision, df_precision.decompression_inmemory, color=colors[2])
    lines!(ax3, df_precision.precision, df_precision.compression_disk, color=colors[3])
    lines!(ax3, df_precision.precision, df_precision.decompression_disk, color=colors[4])

    labels = ["Compression throughput", "Decompression throughput"]
    legend = Legend(fig, [[ln1, ln2], [ln3, ln4]], [labels, labels], ["In-memory", "Disk"], framevisible=false)
    fig[1:3,2] = legend
    display(fig)
    CairoMakie.save("../figs/throughput.pdf", fig)
    return nothing
end

precision_df = precisionTest()
tol_df = tolTest()
rate_df = rateTest()
plot(rate_df, tol_df, precision_df)
