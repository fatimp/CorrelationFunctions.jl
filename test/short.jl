@testset "Test S2 on short predefined sequence" begin
    a = [1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0]

    corr = mean(s2(a, 1; len = 12))
    @test corr ≈ [7, 4, 3, 2, 1, 1, 2, 3, 2, 2, 1, 0] ./ [12 - x for x in 0:11]

    corr = mean(s2(a, 0; len = 12))
    @test corr ≈ [5, 2, 2, 1, 1, 1, 1, 1, 0, 1, 0, 0] ./ [12 - x for x in 0:11]

    corr = mean(s2(a, 1; len = 12, periodic = true))
    @test corr ≈ [7, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4] ./ 12

    corr = mean(s2(a, 0; len = 12, periodic = true))
    @test corr ≈ [5, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2] ./ 12
end

@testset "Test L2 on short predefined sequences" begin
    a = [1, 0, 1, 0, 1, 1, 1]

    corr = mean(l2(a, 1; len = 7))
    @test corr ≈ [5, 2, 1, 0, 0, 0, 0] ./ [7 - x for x in 0:6]

    corr = mean(l2(a, 0; len = 7))
    @test corr ≈ [2, 0, 0, 0, 0, 0, 0] ./ [7 - x for x in 0:6]

    corr = mean(l2(a, 1; len = 7, periodic = true))
    @test corr ≈ [5, 3, 2, 1, 0, 0, 0] ./ 7

    a = [1, 1, 0, 1, 0]
    @test l2(a, 1; len = 5).success ==
        l2(a, 1; len = 5, periodic = true).success

    a = [0, 1, 1, 0, 1]
    @test l2(a, 1; len = 5).success ==
        l2(a, 1; len = 5, periodic = true).success
end

@testset "Test on array with all elements equal to 1" begin
    a = fill(1, (50, 50, 50))
    for p in (false, true)
        for func in (s2, l2)
            @test all(x -> x ≈ 1,
                      mean(func(a, 1; periodic = p, directions = known_directions)))
        end
    end

    for p in (false, true)
        @test all(x -> x ≈ 1,
                  mean(c2(a; periodic = p, directions = known_directions)))
    end

    @test length(pore_size(a, 0).weights) == 0
end
