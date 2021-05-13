diagonals = Directional.diagonals

@testcase "Check diagonal slicer for 2D arrays" begin
    array = Array(reshape(1:40, (5, 8)))
    #5×8 Array{Int64,2}:
    # 1   6  11  16  21  26  31  36
    # 2   7  12  17  22  27  32  37
    # 3   8  13  18  23  28  33  38
    # 4   9  14  19  24  29  34  39
    # 5  10  15  20  25  30  35  40
    diag = collect(diagonals(array, (1, 1)))

    # Expected diagonal slices are
    exp_diag = [[1, 7, 13, 19, 25],
                [2, 8, 14, 20],
                [3, 9, 15],
                [4, 10],
                [5],
                [6, 12, 18, 24, 30],
                [11, 17, 23, 29, 35],
                [16, 22, 28, 34, 40],
                [21, 27, 33, 39],
                [26, 32, 38],
                [31, 37],
                [36]]
    @test diag == exp_diag

    diag = collect(diagonals(array, (-1, 1)))
    # Expected diagonal slices are
    exp_diag = [[1],
                [2, 6],
                [3, 7, 11],
                [4, 8, 12, 16],
                [5, 9, 13, 17, 21],
                [10, 14, 18, 22, 26],
                [15, 19, 23, 27, 31],
                [20, 24, 28, 32, 36],
                [25, 29, 33, 37],
                [30, 34, 38],
                [35, 39],
                [40]]
    @test diag == exp_diag


    array = Array(reshape(1:40, (8, 5)))
    #8×5 Array{Int64,2}:
    # 1   9  17  25  33
    # 2  10  18  26  34
    # 3  11  19  27  35
    # 4  12  20  28  36
    # 5  13  21  29  37
    # 6  14  22  30  38
    # 7  15  23  31  39
    # 8  16  24  32  40
    diag = collect(diagonals(array, (1, 1)))
    exp_diag = [[1, 10, 19, 28, 37],
                [2, 11, 20, 29, 38],
                [3, 12, 21, 30, 39],
                [4, 13, 22, 31, 40],
                [5, 14, 23, 32],
                [6, 15, 24],
                [7, 16],
                [8],
                [9, 18, 27, 36],
                [17, 26, 35],
                [25, 34],
                [33]]
    @test diag == exp_diag

    diag = collect(diagonals(array, (-1, 1)))
    exp_diag = [[1],
                [2, 9],
                [3, 10, 17],
                [4, 11, 18, 25],
                [5, 12, 19, 26, 33],
                [6, 13, 20, 27, 34],
                [7, 14, 21, 28, 35],
                [8, 15, 22, 29, 36],
                [16, 23, 30, 37],
                [24, 31, 38],
                [32, 39],
                [40]]
    @test diag == exp_diag
end

@testcase "Check diagonal slicer for 3D arrays" begin
    # TODO: Maybe something else?
    array = Array(reshape(1:(100^3), (100, 100, 100)))
    # Longest diagonals
    ldiags = (collect(1:10101:1000000),
              collect(100:10099:999901),
              collect(9901:9901:990100),
              collect(990001:-9899:10000))
    directions = (:diag1, :diag2, :diag3, :diag4)

    for (direction, ldiag) in zip(directions, ldiags)
        diags = collect(diagonals(array, Val(direction)))
        flatdiags = reduce(vcat, diags)
        @test length(flatdiags) == length(unique(flatdiags)) == 100^3
        @test ldiag ∈ diags
    end
end
