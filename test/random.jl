# Random array with two phases
rand_array = rand(Float32, (100, 100, 100)) .< 0.3

@testset "Check s²(a, x) = l²(a, x) = 0 for all x where a is an array without phase 2." begin
    for m in (U.Periodic(), U.NonPeriodic())
        for func in (D.s2, D.l2)
            corr = mean_corrfn(func, rand_array, 2; mode = m)
            @test all(corr .== 0)
        end
    end
end

@testset "Check that corr(a, len1, phase) = corr(a, len2, phase)[1:len1] for len2>len1" begin
    for phase in (0, 1)
        for func in (D.s2, D.l2, D.c2, D.surf2, D.surfvoid)
            corr1 = mean_corrfn(func, rand_array, phase; len = 30)
            corr2 = mean_corrfn(func, rand_array, phase; len = 35)
            @test corr1 == corr2[1:30]
        end
    end
end

# Probability of randomly choosing 0 or 1
p = sum(.!rand_array) / length(rand_array) # P{random dot ∈ phase 0}
q = 1 - p                                  # P{random dot ∈ phase 1}
prob = [p, q]

@testset "sᶦ(a,0) = lᶦ(a,0) = P{randomly choosing i}" begin
    for m in (U.Periodic(), U.NonPeriodic())
        for phase in 0:1
            for func in (D.s2, D.l2)
                f = mean_corrfn(func, rand_array, phase; len = 1, mode = m)[1]
                @test relerr(f, prob[phase+1]) < 0.02
            end
        end
    end
end

@testset "sᶦ(a,1) = lᶦ(a,1)" begin
    for m in (U.Periodic(), U.NonPeriodic())
        for phase in 0:1
            s = mean_corrfn(D.s2, rand_array, phase; len = 2, mode = m)[2]
            l = mean_corrfn(D.l2, rand_array, phase; len = 2, mode = m)[2]
            @test s ≈ l
        end
    end
end

@testset "Check some properties of cross_correlation" begin
    for m in (U.Periodic(), U.NonPeriodic())
        cc = D.cross_correlation(rand_array, true, true, U.DirX(); mode = m)
        s2 = D.s2(rand_array, true, U.DirX(); mode = m)
        @test cc == s2

        cc = D.cross_correlation(rand_array, false, false, U.DirX(); mode = m)
        s2 = D.s2(rand_array, false, U.DirX(); mode = m)
        @test cc == s2
    end
end

@testset "Check that l2 is non-increasing" begin
    for m in (U.Periodic(), U.NonPeriodic())
        for phase in 0:1
            for dir in known_directions
                l = D.l2(rand_array, phase, dir; mode = m)
                @test minimum(map(-, l, l[2:end])) >= 0
            end
        end
    end
end

@testset "Check that sⁱ(x) = P{randomly independently choosing i twice} (x > 0)" begin
    for m in (U.Periodic(), U.NonPeriodic())
        for phase in 0:1
            s = mean_corrfn(D.s2, rand_array, phase; mode = m)
            @test relerr(maximum(s[2:end]), prob[phase+1]^2) < 0.04
        end
    end
end

@testset "Check some properties of s3" begin
    ss1, ss2 = U.right_triangles(rand_array, U.PlaneXY())
    for m in (U.Periodic(), U.NonPeriodic())
        s2x = D.s2(rand_array, true, U.DirX(); mode = m)
        s2y = D.s2(rand_array, true, U.DirY(); mode = m)
        s3 = D.s3(rand_array, true, ss1, ss2; mode = m)
        @test all(isapprox.(s3[2:end, 2:end], prob[2]^3; rtol = 0.05))
        @test s3[:,1] ≈ s2x
        @test s3[1,:] ≈ s2y
    end
end

@testset "Check some properties of c3" begin
    ss1, ss2 = U.right_triangles(rand_array, U.PlaneXY())
    for m in (U.Periodic(), U.NonPeriodic())
        c2x = D.c2(rand_array, true, U.DirX(); mode = m)
        c2y = D.c2(rand_array, true, U.DirY(); mode = m)
        c3 = D.c3(rand_array, true, ss1, ss2; mode = m)
        @test c3[:,1] ≈ c2x
        @test c3[1,:] ≈ c2y
    end
end

@testset "Negative shifts" begin
    for _ in 1:100
        p = (rand(0:99), rand(0:99), rand(0:99))
        z = (0, 0, 0)

        ps1 = [z,   z]
        ps2 = [p, .-p]

        for m in (U.Periodic(), U.NonPeriodic())
            s3 = D.s3(rand_array, ps1, ps2; mode = m)
            @test s3[1] == s3[2]
        end
    end
end
