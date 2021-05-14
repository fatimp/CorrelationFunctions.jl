function testreflect(func)
    @testset "$func invariance under mirror transforms" begin
        noise = two_phase_noise_3d()
        for phase in 0:1
            f = mean ∘ func
            # Calculate correlation function on the original image
            corr1 = f(noise, phase)
            # Reflect from yOz plane and calculate correlation function
            corr2 = f(noise[end:-1:1,:,:], phase)
            # Reflect from xOz plane and calculate correlation function
            corr3 = f(noise[:,end:-1:1,:], phase)
            # Reflect from xOy plane and calculate correlation function
            corr4 = f(noise[:,:,end:-1:1], phase)
            @test corr1 ≈ corr2 ≈ corr3 ≈ corr4
        end
    end
end

testreflect(Directional.s2)
testreflect(Directional.l2)
testreflect(Directional.c2)
testreflect(Directional.surfsurf)
# TODO: somehow test histograms returned by pore_size and chord_length

function testsurface(func)
    @testset "Check $(func)⁰(a) = $(func)¹(a) for two phase media" begin
        noise = two_phase_noise_3d()
        @test mean(func(noise, 0; len = 50)) ≈
              mean(func(noise, 1; len = 50))
    end
end

testsurface(Directional.surfsurf)
testsurface(Directional.surfvoid)
