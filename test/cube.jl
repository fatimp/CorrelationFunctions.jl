@testset "Check F_{sss} function for a cube" begin
    cube = zeros(Bool, (300, 300, 300))
    cube[100:200, 100:200, 100:200] .= 1
    @test U.lowfreq_energy_ratio(cube) > 0.97

    ps1 = [(10, 10, 20) for _ in 0:14]
    ps2 = [(20, 20, n)  for n in 0:14]
    pattern = U.ArbitraryPattern(ps1, ps2)
    sss = D.surf3(cube, true, pattern; periodic = true)
    sss = sss * 300^3
    @test all(sss .≈ 2)
end
