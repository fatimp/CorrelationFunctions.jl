@testset "Check distance transform (2D)" begin
    foreach(1:5) do _
        a = BitArray(rand(0:1, (100, 100)))
        @test Utilities.edt(a, Utilities.Plane()) ==
            a |> feature_transform |> distance_transform
    end
end

@testset "Check distance transform (3D)" begin
    foreach(1:5) do _
        a = BitArray(rand(0:1, (100, 100, 100)))
        @test Utilities.edt(a, Utilities.Plane()) ==
            a |> feature_transform |> distance_transform
    end
end

@testset "Check lowfreq_energy_ratio" begin
    a = zeros(Bool, (500, 500))
    a[1:100, 1:100] .= true

    ler  = Utilities.lowfreq_energy_ratio(  a)
    leri = Utilities.lowfreq_energy_ratio(.!a)

    @test ler ≈ leri
    @test ler > Utilities.lowfreq_energy_ratio(rand(Bool, (500, 500)))
end

@testset "Check average_directions" begin
    data = two_phase_noise_3d()
    s2avg  = Directional.s2(data, false; periodic = true) |> mean
    s2avg2 = Map.s2(data, false; periodic = true) |> Map.average_directions
    @test relerr_norm(s2avg, s2avg2) < 0.05
end
