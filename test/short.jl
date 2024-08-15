@testset "Test DIRECTIONAL.S2 on short predefined sequence" begin
    a = [1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0]

    corr = D.s2(a, 1, U.DirX(); len = 12)
    @test corr ≈ [7, 4, 3, 2, 1, 1, 2, 3, 2, 2, 1, 0] ./ [12 - x for x in 0:11]

    corr = D.s2(a, 0, U.DirX(); len = 12)
    @test corr ≈ [5, 2, 2, 1, 1, 1, 1, 1, 0, 1, 0, 0] ./ [12 - x for x in 0:11]

    corr = D.s2(a, 1, U.DirX(); len = 12, mode = U.Periodic())
    @test corr ≈ [7, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4] ./ 12

    corr = D.s2(a, 0, U.DirX(); len = 12, mode = U.Periodic())
    @test corr ≈ [5, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2] ./ 12
end

@testset "Test DIRECTIONAL.L2 on short predefined sequences" begin
    a = [1, 0, 1, 0, 1, 1, 1]

    corr = D.l2(a, 1, U.DirX(); len = 7)
    @test corr ≈ [5, 2, 1, 0, 0, 0, 0] ./ [7 - x for x in 0:6]

    corr = D.l2(a, 0, U.DirX(); len = 7)
    @test corr ≈ [2, 0, 0, 0, 0, 0, 0] ./ [7 - x for x in 0:6]

    corr = D.l2(a, 1, U.DirX(); len = 7, mode = U.Periodic())
    @test corr ≈ [5, 3, 2, 1, 0, 0, 0] ./ 7
end

@testset "Test on array with all elements equal to 1" begin
    a = fill(1, (50, 50, 50))
    for m in (U.Periodic(), U.NonPeriodic())
        for func in (D.s2, D.l2, D.c2)
            for dir in known_directions
                @test all(func(a, 1, dir; mode = m) .≈ 1)
            end
        end
    end

    @test length(D.pore_size(a, 0)) == 0
end
