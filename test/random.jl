# Random array with two phases
rand_array = rand(Float32, (100, 100, 100))
rand_array = map(x -> x < 0.3, rand_array)

@testset "Check s²(a, x) = l²(a, x) = 0 for all x where a is an array without phase 2." begin
    for p in (false, true)
        for func in (D.s2, D.l2)
            corr = mean(func(rand_array, 2; periodic = p, directions = known_directions))
            @test all(x -> x == 0, corr)
        end
    end
end

@testset "Check that corr(a, len1, phase) = corr(a, len2, phase)[1:len1] for len2>len1" begin
    for phase in (0, 1)
        for func in (D.s2, D.l2, D.c2, D.surf2, D.surfvoid)
            corr1 = mean(func(rand_array, phase; len = 30))
            corr2 = mean(func(rand_array, phase; len = 35))
            @test corr1 == corr2[1:30]
        end
    end
end

# Probability of randomly choosing 0 or 1
p = sum(.!rand_array) / length(rand_array) # P{random dot ∈ phase 0}
q = 1 - p                                  # P{random dot ∈ phase 1}
prob = [p, q]

@testset "sᶦ(a,0) = lᶦ(a,0) = P{randomly choosing i}" begin
    for p in (false, true)
        for phase in 0:1
            for func in (D.s2, D.l2)
                f = mean(func(rand_array, phase; len = 1, directions = known_directions))[1]
                @test relerr(f, prob[phase+1]) < 0.02
            end
        end
    end
end

@testset "sᶦ(a,1) = lᶦ(a,1)" begin
    for p in (false, true)
        for phase in 0:1
            s = mean(D.s2(rand_array, phase;
                          len = 2, directions = known_directions))[2]
            l = mean(D.l2(rand_array, phase;
                          len = 2, directions = known_directions))[2]
            @test s ≈ l
        end
    end
end

@testset "Check some properties of cross_correlation" begin
    for periodic in (false, true)
        cc = D.cross_correlation(rand_array, true, true; periodic)
        s2 = D.s2(rand_array, true; periodic)
        @test mean(cc) == mean(s2)

        cc = D.cross_correlation(rand_array, false, false; periodic)
        s2 = D.s2(rand_array, false; periodic)
        @test mean(cc) == mean(s2)
    end
end

@testset "Check that l2 is non-increasing" begin
    for p in (false, true)
        for phase in 0:1
            l = mean(D.l2(rand_array, phase; directions = known_directions))
            @test minimum(map(-, l, l[2:end])) >= 0
        end
    end
end

@testset "Check that sⁱ(x) = P{randomly independently choosing i twice} (x > 0)" begin
    for p in (false, true)
        for phase in 0:1
            s = mean(D.s2(rand_array, phase;
                          directions = known_directions, periodic = p))
            @test relerr(maximum(s[2:end]), prob[phase+1]^2) < 0.04
        end
    end
end

@testset "Check pore size and chord length sum" begin
    for phase in 0:1
        @test sum(D.pore_size(rand_array, phase).weights) ≈ 1
        @test sum(D.chord_length(rand_array, phase).hist.weights) ≈ 1
    end
end

@testset "Check some properties of s3" begin
    for periodic in (false, true)
        s2 = D.s2(rand_array, true; periodic)
        s3 = D.s3(rand_array, true; periodic)
        @test all(isapprox.(s3[D.PlaneXY()][2:end, 2:end], prob[2]^3; rtol = 0.05))
        @test s3[D.PlaneXY()][:,1] ≈ s2[D.DirX()]
        @test s3[D.PlaneXY()][1,:] ≈ s2[D.DirY()]
    end
end

@testset "Check some properties of c3" begin
    for periodic in (false, true)
        c2 = D.c2(rand_array, true; periodic)
        c3 = D.c3(rand_array, true; periodic)
        @test c3[D.PlaneXY()][:,1] ≈ c2[D.DirX()]
        @test c3[D.PlaneXY()][1,:] ≈ c2[D.DirY()]
    end
end
