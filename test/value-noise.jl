function testreflect(func)
    @testset "$func invariance under mirror transforms" begin
        noise = two_phase_noise_3d()
        for phase in 0:1
            for periodic in [false, true]
                # Calculate correlation function on the original image
                corr1 = mean_corrfn(func, noise, phase;
                                    periodic = periodic, directions = axial_directions)
                # Reflect from yOz plane and calculate correlation function
                corr2 = mean_corrfn(func, noise[end:-1:1,:,:], phase;
                                    periodic = periodic, directions = axial_directions)
                # Reflect from xOz plane and calculate correlation function
                corr3 = mean_corrfn(func, noise[:,end:-1:1,:], phase;
                                    periodic = periodic, directions = axial_directions)
                # Reflect from xOy plane and calculate correlation function
                corr4 = mean_corrfn(func, noise[:,:,end:-1:1], phase;
                                    periodic = periodic, directions = axial_directions)
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
        cf1 = mean_corrfn(func, noise, 0; periodic = true, directions = axial_directions)
        cf2 = mean_corrfn(func, noise, 1; periodic = true, directions = axial_directions)
        @test cf1 ≈ cf2
    end
end

testsurface(D.surfvoid)
testsurface(D.surf2)
