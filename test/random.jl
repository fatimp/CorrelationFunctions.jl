const known_directions = [:x,       :y,       :z,
                          :yz_main, :xz_main, :xy_main,
                          :yz_anti, :xz_anti, :xy_anti,
                          :diag1,   :diag2,   :diag3, :diag4]

# Random array with two phases
rand_array = rand(Float32, (50, 51, 52))
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
            corr1 = mfunc(rand_array, 30, phase)
            corr2 = mfunc(rand_array, 35, phase)
            @test corr1 == corr2[1:30]
        end
    end
    corr1 = (mean ∘ c2)(rand_array, 30)
    corr2 = (mean ∘ c2)(rand_array, 35)
    @test corr1 == corr2[1:30]
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
                @test relerr(f, prob[phase+1]) < 0.02
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
            l = mean(l2(rand_array, 40, phase; directions = known_directions))
            @test minimum(map(-, l, l[2:end])) >= 0
        end
    end
end

@testcase "Check that sⁱ(x) = P{randomly independently choosing i twice} (x > 0)" begin
    for p in (false, true)
        for phase in 0:1
            s = mean(s2(rand_array, 40, phase; directions = known_directions, periodic = p))
            @test relerr(maximum(s[2:end]), prob[phase+1]^2) < 0.04
        end
    end
end

@testcase "Check pore size and chord length sum" begin
    for phase in 0:1
        @test sum(pore_size(rand_array, phase).weights) ≈ 1
        @test sum(chord_length(rand_array, phase)[1].weights) ≈ 1
    end
end
