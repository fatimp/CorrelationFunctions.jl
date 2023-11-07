@testset "Check 2D rotations" begin
    q = U.make_rotation(π/2)

    for topology in [U.Torus(), U.Plane()]
        a = rand(Int, (101, 101))
        b = U.rotate_array(a, q, topology)
        @test b == permutedims(a, (2, 1))[:, end:-1:begin]
    end
end

@testset "Check 3D rotations" begin
    q = U.make_rotation(SVector(1, 0, 0), π/2)
    for topology in [U.Torus(), U.Plane()]
        a = rand(Int, (101, 101, 101))
        b = U.rotate_array(a, q, topology)
        @test b == permutedims(a, (1, 3, 2))[:, :, end:-1:begin]
    end

    q = U.make_rotation(SVector(0, 1, 0), π/2)
    for topology in [U.Torus(), U.Plane()]
        a = rand(Int, (101, 101, 101))
        b = U.rotate_array(a, q, topology)
        @test b == permutedims(a, (3, 2, 1))[end:-1:begin, :, :]
    end

    q = U.make_rotation(SVector(0, 0, 1), π/2)
    for topology in [U.Torus(), U.Plane()]
        a = rand(Int, (101, 101, 101))
        b = U.rotate_array(a, q, topology)
        @test b == permutedims(a, (2, 1, 3))[:, end:-1:begin, :]
    end
end

@testset "Check rotation norm" begin
    for _ in 1:100
        rot = U.make_rotation(rand() * 2π)
        @test norm(rot.q) ≈ 1
    end

    for _ in 1:100
        rot = U.make_rotation(SVector(rand(), rand(), rand()), rand() * 2π)
        @test norm(rot.q) ≈ 1
    end
end
