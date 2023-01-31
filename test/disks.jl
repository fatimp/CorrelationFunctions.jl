function gencenters(side, λ)
    n = (λ * side^2) |> Poisson |> rand
    return reduce(hcat, (rand(1:side, 2) for i in 1:n))
end

function gendisks(side, R, λ)
    spheres = zeros(Int8, (side + 2R + 1, side + 2R + 1))
    sphere  = zeros(Int8, (2R + 1, 2R + 1))
    centers = gencenters(side, λ)
    for i in -R:R
        for j in -R:R
            dist = i^2 + j^2
            if dist < R^2
                sphere[j+R+1, i+R+1] = 1
            end
        end
    end
    for center in (centers[:,i] for i in 1:size(centers,2))
        x = center[1]
        y = center[2]
        spheres[x:x + 2R, y:y + 2R] .|= sphere
    end
    return spheres[R+1:end-R-1, R+1:end-R-1]
end

l2_theory(r, R, λ) = exp(-λ*(π*R^2 + 2r*R))
s2_theory(r, R) = (r < 2R) ? 2R^2*acos(r/(2R)) - r/2*sqrt(4R^2 - r^2) : 0
ss_theory(r, R) = (r < 2R) ? 4R^2/(r*sqrt(4R^2-r^2)) : 0
sv_theory(r, R) = (r < 2R) ? 2R*(π - acos(r/(2R))) : 2π*R

pore_size_theory(r, R, λ) = 2λ*π*(r + R)*exp(-λ*π*(r^2 + 2r*R))

@testset "L2 on random overlapping disks generated by Poisson process" begin
    # Area = S^2, radius of a disk = R
    # Poisson parameter = λ
    S = 7000; R = 40; λ = 5e-5; N = 700

    disks = gendisks(S, R, λ)
    calc = D.l2(disks, 0; len = N, periodic = true) |> mean .|> log
    theory = [(log ∘ l2_theory)(r-1, R, λ) for r in 0:N - 1]

    err = relerr.(calc, theory)
    @test maximum(err) < 0.15
end

@testset "S2 for a disk" begin
    S = 500; R = 40

    th(r)  = s2_theory(r, R)
    disk   = draw_ball((S, S), R)
    calc   = D.s2(disk, true; periodic = true) |> mean
    theory = th.(0:length(calc)-1) / S^2

    @test relerr_norm(calc, theory) < 0.01
end

@testset "SS for a disk" begin
    R = 100.2; S = 3000
    boundary = 2R |> floor |> Int

    disk   = draw_ball((S, S), R)
    th(r)  = ss_theory(r, R)
    calc   = D.surfsurf(disk, false; periodic = true) |> mean
    theory = th.(10:boundary-10) / S^2

    @test U.lowfreq_energy_ratio(disk) > 0.97
    @test relerr_norm(calc[10:boundary-10], theory) < 0.08
    @test maximum(calc[boundary+10:end]) < 1e-5
end

@testset "SV for a disk" begin
    S = 3000; R = 100

    th(r)  = sv_theory(r, R)
    disk   = draw_ball((S, S), R)
    calc   = D.surfvoid(disk, false; periodic = true) |> mean
    theory = th.(0:length(calc)-1) / S^2

    @test U.lowfreq_energy_ratio(disk) > 0.97
    @test relerr_norm(calc, theory) < 0.03
end

@testset "Pore size on random overlapping disks generated by Poisson process" begin
    # Area = S^2, radius of a disk = R, Poisson parameter = λ
    S = 7000; R = 50; λ = 5e-5
    disks = gendisks(S, R, λ)

    calc = D.pore_size(disks; nbins = 20)
    edges = calc.edges[1]
    s = step(edges)
    theory = [integrate(x -> pore_size_theory(x, R, λ), n:0.05:n+s)
              for n in 0:s:s*(length(edges) - 2)]

    # Compare cummulative distributions instead of probability
    # densities because density is almost zero for big radii.
    calc_cdf = scan(calc.weights)
    theory_cdf = scan(theory)

    err = relerr.(calc_cdf, theory_cdf)
    @test maximum(err) < 0.1
end

@testset "Chord length for a disk" begin
    S = 500; R = 200

    disk = draw_ball((S, S), R)
    data = D.chord_length(disk, true)
    μ    = π/2 * R
    σ    = sqrt((32 - 3π^2)/12) * R

    @test relerr(data.μ, μ) < 0.02
    @test relerr(data.σ, σ) < 0.02
end
