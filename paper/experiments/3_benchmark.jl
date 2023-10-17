using SequentialCompression
using DataFrames

function compressionthroughput(; rate::Int=0, tol::Real=0, precision::Real=0, inmemory::Bool=true)
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

function compressionthroughput3()
    dummyData = rand(1000,1000,10)
    cp = Vector{Vector{UInt8}}(undef, 10)
    wtime0 = Base.time()
    for it=1:size(dummyData, 3)
        cp[it] = zfp_compress(dummyData[:,:,it])
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

function throughputs(; rate::Int=0, tol::Real=0, precision::Real=0, inmemory::Bool=true)
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
        values[i,:] .= throughputs(rate=rate)
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

function plot(df, title)
    fig = Figure()
    ax = Axis(fig[1,1],
              xlabel="rate (bps)",
              ylabel="throughput (MiB/s)",
              title=title,
              xautolimitmargin=(0,0),
              xticks = (1:4:64)
             )
    ln1 = lines!(ax, df.rate, df.compression_inmemory)
    ln2 = lines!(ax, df.rate, df.decompression_inmemory)
    ln3 = lines!(ax, df.rate, df.compression_disk)
    ln4 = lines!(ax, df.rate, df.decompression_disk)
    labels = ["Compression throughput", "Decompression throughput"]
    legend = Legend(fig, [[ln1, ln2], [ln3, ln4]], [labels, labels], ["In-memory", "Disk"], framevisible=false)
    fig[1,2] = legend
    display(fig)
    return nothing
end

rate_df = rateTest()
plot(rate_df, "Throughput x rate")
