const gcfs = [Map.s2, Map.c2, Map.surfsurf, Map.surfvoid]

const noisegen = [
    () -> rand(Bool, (50)),
    two_phase_noise_2d,
    two_phase_noise_3d,
]

function test_gpu(fn, image)
    for periodic in [true, false]
        cpu_map = fn(image,          true; periodic = periodic)
        gpu_map = fn(CuArray(image), true; periodic = periodic)
        @test CuArray(cpu_map) ≈ gpu_map
    end
end

for noise_fn in noisegen
    image = noise_fn()
    for fn in gcfs
        test_gpu(fn, image)
    end
end