const noisegen = [
    () -> rand(Bool, (50)),
    two_phase_noise_2d,
    two_phase_noise_3d,
]

function test_gpu(fn, image)
    for periodic in [true, false]
        # FIXME: only periodic boundary conditions for filters are
        # supported now.
        local cpu_map, gpu_map
        if fn ∈ [M.surfsurf, M.surfvoid]
            flt = U.EdgeFilter(U.BCPeriodic(), U.Kernel5x5())
            cpu_map = fn(image,          true;
                         periodic = periodic, filter = flt)
            gpu_map = fn(CuArray(image), true;
                         periodic = periodic, filter = flt)
        else
            cpu_map = fn(image,          true; periodic = periodic)
            gpu_map = fn(CuArray(image), true; periodic = periodic)
        end
        @test CuArray(cpu_map) ≈ gpu_map
    end
end

for noise_fn in noisegen
    image = noise_fn()
    for fn in cfs[ndims(image)]
        test_gpu(M.eval(fn), image)
    end
end
