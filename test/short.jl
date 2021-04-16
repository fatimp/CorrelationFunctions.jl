@testset "Test S2 on short predefined sequence" begin
    a = [1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0]

    corr = s2(a, 1; len = 12)
    @test corr.success[:x] == [7, 4, 3, 2, 1, 1, 2, 3, 2, 2, 1, 0]
    @test corr.total[:x] == [12 - x for x in 0:11]

    corr = s2(a, 0; len = 12)
    @test corr.success[:x] == [5, 2, 2, 1, 1, 1, 1, 1, 0, 1, 0, 0]
    @test corr.total[:x] == [12 - x for x in 0:11]

    corr = s2(a, 1; len = 12, periodic = true)
    @test corr.success[:x] == [7, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4]
    @test corr.total[:x] == [12 for x in 0:11]

    corr = s2(a, 0; len = 12, periodic = true)
    @test corr.success[:x] == [5, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2]
    @test corr.total[:x] == [12 for x in 0:11]
end

@testset "Test L2 on short predefined sequences" begin
    a = [1, 0, 1, 0, 1, 1, 1]

    corr = l2(a, 1; len = 7)
    @test corr.success[:x] == [5, 2, 1, 0, 0, 0, 0]
    @test corr.total[:x] == [7 - x for x in 0:6]

    corr = l2(a, 0; len = 7)
    @test corr.success[:x] == [2, 0, 0, 0, 0, 0, 0]
    @test corr.total[:x] == [7 - x for x in 0:6]

    corr = l2(a, 1; len = 7, periodic = true)
    @test corr.success[:x] == [5, 3, 2, 1, 0, 0, 0]
    @test corr.total[:x] == [7 for x in 0:6]

    a = [1, 1, 0, 1, 0]
    @test l2(a, 1; len = 5).success ==
          l2(a, 1; len = 5, periodic = true).success

    a = [0, 1, 1, 0, 1]
    @test l2(a, 1; len = 5).success ==
          l2(a, 1; len = 5, periodic = true).success
end
