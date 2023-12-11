const USE_GPU = false  # Use GPU? If this is set false, then no GPU needs to be available
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
using SequentialCompression
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using Printf, Statistics
using JSON

@parallel function compute_V!(Vx::Data.Array, Vy::Data.Array, P::Data.Array, dt::Data.Number, ρ::Data.Number, dx::Data.Number, dy::Data.Number)
    @inn(Vx) = @inn(Vx) - dt/ρ*@d_xi(P)/dx
    @inn(Vy) = @inn(Vy) - dt/ρ*@d_yi(P)/dy
    return
end

@parallel function compute_P!(P::Data.Array, Vx::Data.Array, Vy::Data.Array, dt::Data.Number, k::Data.Number, dx::Data.Number, dy::Data.Number)
    @all(P) = @all(P) - dt*k*(@d_xa(Vx)/dx + @d_ya(Vy)/dy)
    return
end

##################################################
@views function acoustic2D()
    # Physics
    lx, ly    = 50.0, 50.0  # domain extends
    k         = 1.0         # bulk modulus
    ρ         = 1.0         # density
    t         = 0.0         # physical time
    # Numerics
    nx, ny    = 255, 255    # numerical grid resolution; should be a mulitple of 32-1 for optimal GPU perf
    nt        = 1000        # number of timesteps
    nout      = 10          # plotting frequency
    # Derived numerics
    dx, dy    = lx/(nx-1), ly/(ny-1) # cell sizes
    dx, dy = map(d->round(d, digits=2), (dx, dy))
    # Array allocations
    P         = @zeros(nx  ,ny  )
    Vx        = @zeros(nx+1,ny  )
    Vy        = @zeros(nx  ,ny+1)
    # Initial conditions
    P        .= Data.Array([exp(-((ix-1)*dx-0.5*lx)^2 -((iy-1)*dy-0.5*ly)^2) for ix=1:size(P,1), iy=1:size(P,2)])
    dt        = min(dx,dy)/sqrt(k/ρ)/4.1

    snapshots = SeqCompressor(Float64, nx, ny, inmemory=false)

    # Time loop
    for it = 1:nt
        if (it==11)  global wtime0 = Base.time()  end
        @parallel compute_V!(Vx, Vy, P, dt, ρ, dx, dy)
        @parallel compute_P!(P, Vx, Vy, dt, k, dx, dy)
        append!(snapshots, P)
        t = t + dt
    end
    # Performance
    wtime    = Base.time()-wtime0
    A_eff    = (3*2)/1e9*nx*ny*sizeof(Data.Number)  # Effective main memory access per iteration [GB] (Lower bound of required memory access: H and dHdτ have to be read and written (dHdτ for damping): 4 whole-array memaccess; B has to be read: 1 whole-array memaccess)
    wtime_it = wtime/(nt-10)                        # Execution time per iteration [s]
    T_eff    = A_eff/wtime_it                       # Effective memory throughput [GB/s]
    @printf("Total steps=%d, time=%1.3e sec (@ T_eff = %1.2f GB/s) \n", nt, wtime, round(T_eff, sigdigits=2))

    SequentialCompression.save("wavefield.szfp", snapshots)

    param = Dict("dt" => dt, "dx" => dx, "dy" => dy)
    open("geometry.json", "w") do file
        JSON.print(file, param)
    end

    return
end

acoustic2D()
