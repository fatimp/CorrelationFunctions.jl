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
        l = (mean ∘ Directional.l2)(cb, phase)
        @test l[1] ≈ 1/2
        # Would be 1/4 exactly on infinite checkboard
        @test isapprox(l[2], 1/4; atol = 1e-2)
        @test all(x -> x == 0, l[3:end])
    end
end

# C2 behaves exactly as L2 with segmentation algorithm chosen which
# does not join segments in diagonal directions.
@testset "C2" begin
    c = (mean ∘ Directional.c2)(cb, 1)
    @test c[1] ≈ 1/2
    # Would be 1/4 exactly on infinite checkboard
    @test isapprox(c[2], 1/4; atol = 10^-2)
    @test all(x -> x == 0, c[3:end])
end

@testset "S2" begin
    for phase in 0:1
        s = mean(Directional.s2(cb, phase; periodic = true))
        @test all(x -> x ≈ 1/4, s[2:2:end])
        @test all(x -> isapprox(x, 0; atol = 1e-2), s[3:4:end])
        @test all(x -> x ≈ 1/2, s[1:4:end])
    end
end
