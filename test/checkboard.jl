"""
    gen_checkboard(side)

Generate cube array of ones and zeros with checkboard-like
pattern. The size of checkboard cell is 2.
"""
function gen_checkboard(side)
    checkboard = Array{Int, 3}(undef, side, side, side)
    for i in 1:side
        for j in 1:side
            for k in 1:side
                x = ((i-1) & 2 != 0) ? 1 : 0
                y = ((j-1) & 2 != 0) ? 1 : 0
                z = ((k-1) & 2 != 0) ? 1 : 0
                checkboard[k,j,i] = x ⊻ y ⊻ z
            end
        end
    end
    return checkboard
end

cb = gen_checkboard(60)

@testset "L2" begin
    for phase in 0:1
        l = mean_corrfn(D.l2, cb, phase;
                        periodic = true, directions = axial_directions)
        @test l[1] ≈ 1/2
        @test l[2] ≈ 1/4
        @test all(l[3:end] .== 0)
    end
end

# C2 behaves exactly as L2 with segmentation algorithm chosen which
# does not join segments in diagonal directions.
@testset "C2" begin
    for phase in 0:1
        c = mean_corrfn(D.c2, cb, phase;
                        periodic = true, directions = axial_directions)
        @test c[1] ≈ 1/2
        @test c[2] ≈ 1/4
        @test all(c[3:end] .== 0)
    end
end

@testset "S2" begin
    for phase in 0:1
        s = mean_corrfn(D.s2, cb, phase;
                        periodic = true, directions = axial_directions)
        @test all(s[2:2:end] .≈ 1/4)
        @test all(x -> isapprox(x, 0; atol = 1e-2), s[3:4:end])
        @test all(s[1:4:end] .≈ 1/2)
    end
end
