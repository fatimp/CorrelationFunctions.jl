macro testreflect(func)
    quote
        @testset $("$func invariance under mirror transforms") begin
            noise = two_phase_noise()
            for phase in 0:1
                f = mean ∘ $func
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
end

@testreflect s2
@testreflect l2
@testreflect c2
@testreflect surfsurf
# TODO: somehow test histograms returned by pore_size and chord_length

macro testsurface(func)
    quote
        @testset $("Check $(func)⁰(a) = $(func)¹(a) for two phase media") begin
            noise = two_phase_noise()
            @test mean($func(noise, 0; len = 50)) ≈
                  mean($func(noise, 1; len = 50))
        end
    end
end

@testsurface surfsurf
@testsurface surfvoid
