function testreflect(func)
    @testset "$func invariance under mirror transforms" begin
        noise = two_phase_noise_3d()
        for phase in 0:1
            for periodic in [false, true]
                f = mean ∘ func
                # Calculate correlation function on the original image
                corr1 = f(noise, phase; periodic = periodic)
                # Reflect from yOz plane and calculate correlation function
                corr2 = f(noise[end:-1:1,:,:], phase; periodic = periodic)
                # Reflect from xOz plane and calculate correlation function
                corr3 = f(noise[:,end:-1:1,:], phase; periodic = periodic)
                # Reflect from xOy plane and calculate correlation function
                corr4 = f(noise[:,:,end:-1:1], phase; periodic = periodic)
                @test corr1 ≈ corr2 ≈ corr3 ≈ corr4
            end
        end
    end
end

ssfp(array, phase; periodic) =
    Directional.surfsurf(array, phase;
                         periodic, filter = Utilities.EdgeFilter(Utilities.edge_5x5,
                                                                 Utilities.BCPeriodic()))
ssfz(array, phase; periodic) =
    Directional.surfsurf(array, phase;
                         periodic, filter = Utilities.EdgeFilter(Utilities.edge_5x5,
                                                                 Utilities.BCReflect()))
svfp(array, phase; periodic) =
    Directional.surfvoid(array, phase;
                         periodic, filter = Utilities.EdgeFilter(Utilities.edge_5x5,
                                                                 Utilities.BCPeriodic()))
svfz(array, phase; periodic) =
    Directional.surfvoid(array, phase;
                         periodic, filter = Utilities.EdgeFilter(Utilities.edge_5x5,
                                                                 Utilities.BCReflect()))

testreflect(Directional.s2)
testreflect(Directional.l2)
testreflect(Directional.c2)
testreflect(ssfp)
testreflect(ssfz)
# TODO: somehow test histograms returned by pore_size and chord_length

function testsurface(func)
    @testset "Check $(func)⁰(a) = $(func)¹(a) for two phase media" begin
        noise = two_phase_noise_3d()
        @test mean(func(noise, 0; periodic = true)) ≈ mean(func(noise, 1; periodic = true))
    end
end

testsurface(ssfp)
testsurface(ssfz)
testsurface(svfp)
testsurface(svfz)
