const known_directions = [:x,       :y,       :z,
                          :yz_main, :xz_main, :xy_main,
                          :yz_anti, :xz_anti, :xy_anti,
                          :diag1,   :diag2,   :diag3, :diag4]

# Random array with two phases
rand_array = rand(Float32, (70, 80, 90))
rand_array = map(x -> (x<0.3) ? 0 : 1, rand_array)

@testcase "Check s²(a, x) = l²(a, x) = 0 for all x where a is an array without phase 2." begin
    for p in (false, true)
        for func in (s2, l2)
            corr = mean(func(rand_array, 40, 2; periodic = p, directions = known_directions))
            @test all(x -> x == 0, corr)
        end
    end
end

@testcase "Check that corr(a, len1, phase) = corr(a, len2, phase)[1:len1] for len2>len1" begin
    for phase in (0, 1)
        for func in (s2, l2, surfsurf)
            mfunc = mean ∘ func
            corr1 = mfunc(rand_array, 40, phase)
            corr2 = mfunc(rand_array, 70, phase)
            @test corr1 == corr2[1:40]
        end
    end
    corr1 = (mean ∘ c2)(rand_array, 40)
    corr2 = (mean ∘ c2)(rand_array, 70)
    @test corr1 == corr2[1:40]
end

# Probability of randomly choosing 0 or 1
p = count(x -> x == 0, rand_array) / length(rand_array) # P{random dot ∈ phase 0}
q = 1 - p                                               # P{random dot ∈ phase 1}
prob = [p, q]

@testcase "sᶦ(a,0) = lᶦ(a,0) = P{randomly choosing i}" begin
    for p in (false, true)
        for phase in 0:1
            for func in (s2, l2)
                f = mean(func(rand_array, 1, phase, directions = known_directions))[1]
                @test abs(f - prob[phase+1]) < 1e-3
            end
        end
    end
end

@testcase "sᶦ(a,1) = lᶦ(a,1)" begin
    for p in (false, true)
        for phase in 0:1
            s = mean(s2(rand_array, 2, phase, directions = known_directions))[2]
            l = mean(l2(rand_array, 2, phase, directions = known_directions))[2]
            @test s ≈ l
        end
    end
end

@testcase "Check that l2 is non-increasing" begin
    for p in (false, true)
        for phase in 0:1
            l = mean(l2(rand_array, 50, phase; directions = known_directions))
            @test all(x -> x >= 0, map(-, l, l[2:end]))
        end
    end
end

@testcase "Check that sⁱ(x) = P{randomly independently choosing i twice} (x > 0)" begin
    for phase in 0:1
        s = mean(s2(rand_array, 50, phase; directions = known_directions))
        @test all(x -> abs(x-prob[phase+1]^2) < 1e-3, s[2:end])
    end

    # FIXME: check periodic case
    for phase in 0:1
        s = mean(s2(rand_array, 50, phase; directions = known_directions, periodic = true))
        @test_broken all(x -> abs(x-prob[phase+1]^2) < 1e-3, s[2:end])
    end
end
