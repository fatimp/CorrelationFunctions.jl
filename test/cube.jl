@testset "Check F_{sss} function for a cube" begin
    cube = zeros(Bool, (300, 300, 300))
    cube[100:200, 100:200, 100:200] .= 1
    @test U.lowfreq_energy_ratio(cube) > 0.97

    sss = D.surf3(cube, true, [(10, 10, 20)], [(20, 20, n) for n in 0:14]; periodic = true)
    sss = sss * 300^3
    @test all(sss .≈ 2)
end
