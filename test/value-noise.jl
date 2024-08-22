function testreflect(func)
    @testset "$func invariance under mirror transforms" begin
        noise = two_phase_noise_3d()
        for phase in 0:1
            for m in (U.Periodic(), U.NonPeriodic())
                # Calculate correlation function on the original image
                corr1 = mean_corrfn(func, noise, phase;
                                    mode = m, directions = axial_directions)
                # Reflect from yOz plane and calculate correlation function
                corr2 = mean_corrfn(func, noise[end:-1:1,:,:], phase;
                                    mode = m, directions = axial_directions)
                # Reflect from xOz plane and calculate correlation function
                corr3 = mean_corrfn(func, noise[:,end:-1:1,:], phase;
                                    mode = m, directions = axial_directions)
                # Reflect from xOy plane and calculate correlation function
                corr4 = mean_corrfn(func, noise[:,:,end:-1:1], phase;
                                    mode = m, directions = axial_directions)
                @test corr1 ≈ corr2 ≈ corr3 ≈ corr4
            end
        end
    end
end

testreflect(D.s2)
testreflect(D.l2)
testreflect(D.c2)
testreflect(D.surf2)
# TODO: somehow test histograms returned by pore_size and chord_length

function testsurface(func)
    @testset "Check $(func)⁰(a) = $(func)¹(a) for two phase media" begin
        noise = two_phase_noise_3d()
        cf1 = mean_corrfn(func, noise, 0; mode = U.Periodic(), directions = axial_directions)
        cf2 = mean_corrfn(func, noise, 1; mode = U.Periodic(), directions = axial_directions)
        @test cf1 ≈ cf2
    end
end

testsurface(D.surfvoid)
testsurface(D.surf2)

@testset "Check S₂ calculation with a mask (Map)" begin
    noise = two_phase_noise_3d()
    big = rand(Bool, (100, 100, 100))
    big[1:50, 1:50, 1:50] = noise
    mask = zeros(Bool, (100, 100, 100))
    mask[1:50, 1:50, 1:50] .= 1

    for phase in [false, true]
        cf1 = M.s2(big, phase; mode = U.Mask(mask))
        cf2 = M.s2(noise, phase; mode = U.NonPeriodic())
        @test cf1[1:50, 1:50, 1:50] ≈ cf2[1:50, 1:50, 1:50]
    end
end

@testset "Check S₂ calculation with a mask (Directional)" begin
    noise = two_phase_noise_3d()
    big = rand(Bool, (100, 100, 100))
    big[1:50, 1:50, 1:50] = noise
    mask = zeros(Bool, (100, 100, 100))
    mask[1:50, 1:50, 1:50] .= 1

    for phase in [false, true]
        cf1 = D.s2(big,   phase, U.DirX(); len = 50, mode = U.Mask(mask))
        cf2 = D.s2(noise, phase, U.DirX(); len = 50, mode = U.NonPeriodic())
        @test cf1 ≈ cf2
    end
end

@testset "Check CC calculation with a mask (Map)" begin
    noise = two_phase_noise_3d()
    big = rand(Bool, (100, 100, 100))
    big[1:50, 1:50, 1:50] = noise
    mask = zeros(Bool, (100, 100, 100))
    mask[1:50, 1:50, 1:50] .= 1

    cf1 = M.cross_correlation(big, false, true; mode = U.Mask(mask))
    cf2 = M.cross_correlation(noise, false, true; mode = U.NonPeriodic())
    @test cf1[1:50, 1:50, 1:50] ≈ cf2[1:50, 1:50, 1:50]
end

@testset "Check CC calculation with a mask (Directional)" begin
    noise = two_phase_noise_3d()
    big = rand(Bool, (100, 100, 100))
    big[1:50, 1:50, 1:50] = noise
    mask = zeros(Bool, (100, 100, 100))
    mask[1:50, 1:50, 1:50] .= 1

    cf1 = D.cross_correlation(big,   false, true, U.DirX(); len = 50, mode = U.Mask(mask))
    cf2 = D.cross_correlation(noise, false, true, U.DirX(); len = 50, mode = U.NonPeriodic())
    @test cf1 ≈ cf2
end

@testset "Check F_{ss} calculation with a mask (Map)" begin
    noise = two_phase_noise_3d()
    # Add zero padding because of different boundary conditions in extract_edges()
    padded = zeros(Bool, (60, 60, 60))
    padded[6:55, 6:55, 6:55] = noise
    big = rand(Bool, (100, 100, 100))
    big[1:60, 1:60, 1:60] = padded
    mask = zeros(Bool, (100, 100, 100))
    mask[1:60, 1:60, 1:60] .= 1

    for phase in [false, true]
        cf1 = M.surf2(big, phase; mode = U.Mask(mask))
        cf2 = M.surf2(padded, phase; mode = U.NonPeriodic())
        @test cf1[1:60, 1:60, 1:60] ≈ cf2[1:60, 1:60, 1:60]
    end
end

#=
@testset "Check F_{ss} calculation with a mask (Directional)" begin
    noise = two_phase_noise_3d()
    # Add zero padding because of different boundary conditions in extract_edges()
    padded = zeros(Bool, (60, 60, 60))
    padded[6:55, 6:55, 6:55] = noise
    big = rand(Bool, (100, 100, 100))
    big[1:60, 1:60, 1:60] = padded
    mask = zeros(Bool, (100, 100, 100))
    mask[1:60, 1:60, 1:60] .= 1

    for phase in [false, true]
        cf1 = D.surf2(big,    phase, U.DirX(); len = 60, mode = U.Mask(mask))
        cf2 = D.surf2(padded, phase, U.DirY(); len = 60, mode = U.NonPeriodic())
        @test isapprox(cf1, cf2; atol = 1e-3)
        #@test cf1 ≈ cf2
    end
end
=#

@testset "Check F_{sv} calculation with a mask (Map)" begin
    noise = two_phase_noise_3d()
    # Add zero padding because of different boundary conditions in extract_edges()
    padded = zeros(Bool, (60, 60, 60))
    padded[6:55, 6:55, 6:55] = noise
    big = rand(Bool, (100, 100, 100))
    big[1:60, 1:60, 1:60] = padded
    mask = zeros(Bool, (100, 100, 100))
    mask[1:60, 1:60, 1:60] .= 1

    for phase in [false, true]
        cf1 = M.surfvoid(big, phase; mode = U.Mask(mask))
        cf2 = M.surfvoid(padded, phase; mode = U.NonPeriodic())
        @test cf1[1:60, 1:60, 1:60] ≈ cf2[1:60, 1:60, 1:60]
    end
end

@testset "Check F_{sv} calculation with a mask (Directional)" begin
    noise = two_phase_noise_3d()
    # Add zero padding because of different boundary conditions in extract_edges()
    padded = zeros(Bool, (60, 60, 60))
    padded[6:55, 6:55, 6:55] = noise
    big = rand(Bool, (100, 100, 100))
    big[1:60, 1:60, 1:60] = padded
    mask = zeros(Bool, (100, 100, 100))
    mask[1:60, 1:60, 1:60] .= 1

    for phase in [false, true]
        cf1 = D.surfvoid(big,    phase, U.DirX(); len = 60, mode = U.Mask(mask))
        cf2 = D.surfvoid(padded, phase, U.DirX(); len = 60, mode = U.NonPeriodic())
        @test cf1 ≈ cf2
    end
end

@testset "Check C₂ calculation with a mask (Map)" begin
    noise = two_phase_noise_3d()[1:20, 1:20, 1:20] # Slow otherwise
    big = rand(Bool, (50, 50, 50))
    big[1:20, 1:20, 1:20] = noise
    mask = zeros(Bool, (50, 50, 50))
    mask[1:20, 1:20, 1:20] .= 1

    for phase in [false, true]
        cf1 = M.c2(big, phase; mode = U.Mask(mask))
        cf2 = M.c2(noise, phase; mode = U.NonPeriodic())
        @test cf1[1:20, 1:20, 1:20] ≈ cf2[1:20, 1:20, 1:20]
    end
end

@testset "Check C₂ calculation with a mask (Directional)" begin
    noise = two_phase_noise_3d()[1:20, 1:20, 1:20] # Slow otherwise
    big = rand(Bool, (50, 50, 50))
    big[1:20, 1:20, 1:20] = noise
    mask = zeros(Bool, (50, 50, 50))
    mask[1:20, 1:20, 1:20] .= 1

    for phase in [false, true]
        cf1 = D.c2(big,   phase, U.DirX(); len = 20, mode = U.Mask(mask))
        cf2 = D.c2(noise, phase, U.DirX(); len = 20, mode = U.NonPeriodic())
        @test cf1 ≈ cf2
    end
end
