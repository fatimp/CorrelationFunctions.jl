testdist(x, y, :: U.NonPeriodic) = norm(x - y)

function testdist(x, y, :: U.Periodic)
    c = map(x, y) do i, j
        d1 = abs(i - j)
        d2 = 100 - d1
        d = min(d1, d2)
    end

    return norm(c)
end

function check_edt(dim, topology)
    toarray(x) = Int[Tuple(x)...]
    foreach(1:5) do _
        dims = (100 for _ in 1:dim)
        a = falses(dims...)
        idx1 = CartesianIndex((rand(1:100) for _ in 1:dim)...)
        idx2 = CartesianIndex((rand(1:100) for _ in 1:dim)...)
        a[idx1] = 1
        a[idx2] = 1

        edt = map(CartesianIndices(a)) do idx
            d1 = testdist(toarray(idx), toarray(idx1), topology)
            d2 = testdist(toarray(idx), toarray(idx2), topology)
            min(d1, d2)
        end

        @test U.edt(a, topology) ≈ edt
    end
end

@testset "Check distance transform (2D)" begin
    check_edt(2, U.Periodic())
    check_edt(2, U.NonPeriodic())
end

@testset "Check distance transform (3D)" begin
    check_edt(3, U.Periodic())
    check_edt(3, U.NonPeriodic())
end

@testset "Check lowfreq_energy_ratio" begin
    a = zeros(Bool, (500, 500))
    a[1:100, 1:100] .= true

    ler  = U.lowfreq_energy_ratio(  a)
    leri = U.lowfreq_energy_ratio(.!a)

    @test ler ≈ leri
    @test ler > U.lowfreq_energy_ratio(rand(Bool, (500, 500)))
end

@testset "Check average_directions" begin
    data = two_phase_noise_3d()
    s2avg  = mean_corrfn(D.s2, data, false;
                         mode = U.Periodic(), directions = known_directions)
    s2avg2 = M.s2(data, false; mode = U.Periodic()) |> M.average_directions
    @test relerr_norm(s2avg, s2avg2) < 0.05
end
