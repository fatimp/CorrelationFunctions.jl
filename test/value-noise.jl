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

@testset "Check S₂ calculation with trivial mask" begin
    noise = two_phase_noise_3d()
    mask = ones(Bool, size(noise))

    for phase in [false, true]
        @test M.s2(noise, phase; mode = U.Mask(mask)) ≈
            M.s2(noise, phase; mode = U.NonPeriodic())
    end
end

@testset "Check CC calculation with trivial mask" begin
    noise = two_phase_noise_3d()
    mask = ones(Bool, size(noise))

    @test M.cross_correlation(noise, true, false; mode = U.Mask(mask)) ≈
        M.cross_correlation(noise, true, false; mode = U.NonPeriodic())
    @test M.cross_correlation(noise, false, true; mode = U.Mask(mask)) ≈
        M.cross_correlation(noise, false, true; mode = U.NonPeriodic())
end

@testset "Check C₂ calculation with trivial mask" begin
    noise = two_phase_noise_3d()
    mask = ones(Bool, size(noise))

    for phase in [true, false]
        @test M.c2(noise, phase; mode = U.Mask(mask)) ≈
            M.c2(noise, phase; mode = U.NonPeriodic())
    end
end
