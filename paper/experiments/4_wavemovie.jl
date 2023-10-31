using GLMakie
using SequentialCompression


function makemovie(A::SequentialCompression.AbstractCompArraySeq, videoFilePath::String)
    fig = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98))
    ax = Axis(fig[1,1])

    iter = Observable(1)
    img = @lift(image!(ax, A[$iter]))

    framerate = 30

    iterations = 1:size(A)[3]

    display(fig)

    if !occursin(".mkv", videoFilePath)
        videoFilePath *= ".mkv"
    end

    record(fig, videoFilePath, iterations;
           framerate=framerate) do it
        iter[] = it
    end

    return nothing
end

A = load("./wavefield.szfp")
makemovie(A, "wavemovie")
