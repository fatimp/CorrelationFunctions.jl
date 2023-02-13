const noisegen = [
    () -> rand(Bool, (50)),
    two_phase_noise_2d,
    two_phase_noise_3d,
]

function test_surf(fn, image)
    for kernel in [U.ConvKernel(5), U.ErosionKernel(5)]
        cpu_map = fn(image, true;          periodic = true, filter = kernel)
        gpu_map = fn(CuArray(image), true; periodic = true, filter = kernel)
        @test CuArray(cpu_map) ≈ gpu_map
    end
end

function test_prob(fn, image)
    for periodic in [true, false]
        cpu_map = fn(image,          true; periodic)
        gpu_map = fn(CuArray(image), true; periodic)
        @test CuArray(cpu_map) ≈ gpu_map
    end
end

for noise_fn in noisegen
    image = noise_fn()
    for fn in cfs[ndims(image)]
        func = M.eval(fn)
        if func ∈ [M.surfsurf, M.surfvoid]
            test_surf(func, image)
        else
            test_prob(func, image)
        end
    end
end
