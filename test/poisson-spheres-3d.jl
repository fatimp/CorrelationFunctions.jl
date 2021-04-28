function gencenters(side, λ)
    n = pois_rand(λ * side^3)
    return reduce(hcat, (rand(1:side, 3) for i in 1:n))
end

function genballs(side, R, λ)
    sphere_side = 2R + 1
    padded_side = side + sphere_side
    spheres = zeros(Int8, (padded_side, padded_side, padded_side))
    sphere = zeros(Int8, (sphere_side, sphere_side, sphere_side))
    centers = gencenters(side, λ)
    for i in -R:R
        for j in -R:R
            for k in -R:R
                dist = i^2 + j^2 + k^2
                if dist < R^2
                    sphere[k+R+1, j+R+1, i+R+1] = 1
                end
            end
        end
    end

    for center in (centers[:,i] for i in 1:size(centers,2))
        x = center[1]
        y = center[2]
        z = center[3]
        spheres[x:x + 2R, y:y + 2R, z:z + 2R] .|= sphere
    end
    return spheres[R+1:end-R-1, R+1:end-R-1, R+1:end-R-1]
end

heaviside(x) = max(sign(x), 0)

function s2_theory(r, R, λ)
    tmp = r/R
    tmp2 = (r > 2R) ? 2 : 1 + 3/4*tmp - 1/16*tmp^3
    η = λ * 4/3 * π * R^3
    return exp(-η*tmp2)
end

function ss_theory(r, R, λ)
    η = λ * 4/3 * π * R^3
    tmp = 3η/R
    tmp2 = (tmp * (1 - (0.5 - r/(4R))*heaviside(2R-r)))^2 + tmp/(2r)*heaviside(2R-r)
    return s2_theory(r, R, λ) * tmp2
end

function sv_theory(r, R, λ)
    η = λ * 4/3 * π * R^3
    tmp = 3η/R
    return s2_theory(r, R, λ) * tmp * (1 - (0.5 - r/(4R))*heaviside(2R-r))
end

function l2_theory(r, R, λ)
    η = λ * 4/3 * π * R^3
    p = exp(-η)
    return p^(1 + 3r/(4R))
end

function pore_size_theory(r, R, λ)
    s(x) = 4π * x^2
    v(x) = 4/3 * π * x^3
    η  = v(R) * λ
    p  = exp(-η)
    s1 = s(r + R)
    v1 = v(r + R)
    return λ * s1 * exp(-λ * v1) / p
end

function chord_length_theory(r, R, λ)
    η = λ * 4/3 * π * R^3
    p = exp(-η)
    return 3η/(4R)*p^(3r/(4R))
end

mean_chord_length(R, λ) = 1/(λ * π *R^2)

@testset "S2 on random overlapping balls generated by Poisson process" begin
    # Volume = S^3, radius of a ball = R
    # Poisson parameter = λ
    S = 500; R = 20; λ = 3e-5
    balls = genballs(S, R, λ)

    # Calculate in one direction for speed
    calc = mean(s2(balls, 0; directions = [:x]))
    theory = [s2_theory(r - 1, R, λ) for r in 0:length(calc) - 1]
    err = relerr.(calc, theory)

    @test maximum(err) < 0.25
end

@testset "SS on random overlapping balls generated by Poisson process" begin
    # Volume = S^3, radius of a ball = R
    # Poisson parameter = λ
    S = 500; R = 20; λ = 3e-5;
    balls = genballs(S, R, λ)

    # Calculate in one direction for speed
    calc = mean(surfsurf(balls, 0; directions = [:x]))
    theory = [ss_theory(r - 1, R, λ) for r in 0:length(calc) - 1]
    err = relerr.(calc, theory)

    # Surface-surface function with r < 2R has a discontinuity and two
    # points where it goes to +Inf. Do not test it.
    @test maximum(err[2R + 10:end]) < 0.25
end

@testset "SV on random overlapping balls generated by Poisson process" begin
    # Volume = S^3, radius of a ball = R
    # Poisson parameter = λ
    S = 500; R = 20; λ = 3e-5;
    balls = genballs(S, R, λ)

    # Calculate in one direction for speed
    calc = mean(surfvoid(balls, 0; directions = [:x]))
    theory = [sv_theory(r - 1, R, λ) for r in 0:length(calc) - 1]
    err = relerr.(calc, theory)

    @test maximum(err) < 0.25
end

@testset "L2 on random overlapping balls generated by Poisson process" begin
    # Volume = S^3, radius of a ball = R
    # Poisson parameter = λ
    # Compare l2(x) for x ∈ {0,1,...,N-1}
    S = 500; R = 20; λ = 3e-5; N = 100
    balls = genballs(S, R, λ)

    calc = log.(mean(l2(balls, 0; len = N)))
    theory = [log(l2_theory(r - 1, R, λ)) for r in 0:N - 1]
    err = relerr.(calc, theory)

    @test maximum(err) < 0.25
end

@testset "Pore size on random overlapping balls generated by Poisson process" begin
    # Volume = S^3, radius of a ball = R
    # Poisson parameter = λ
    S = 400; R = 30; λ = 1e-5
    balls = genballs(S, R, λ)

    calc = pore_size(balls; nbins = 20)
    edges = calc.edges[1]
    s = step(edges)
    theory = [integrate(x -> pore_size_theory(x, R, λ), n:0.05:n+s)
              for n in 0:s:s*(length(edges) - 2)]

    # Compare cummulative distributions instead of probability
    # densities because density is almost zero for big radii.
    calc_cdf = scan(calc.weights)
    theory_cdf = scan(theory)

    err = relerr.(calc_cdf, theory_cdf)
    @test maximum(err) < 0.25
end

@testset "Chord length on random overlapping balls generated by Poisson process" begin
    # Volume = S^3, radius of a ball = R
    # Poisson parameter = λ
    S = 400; R = 30; λ = 1e-5
    balls = genballs(S, R, λ)

    calc, mc = chord_length(balls, 0; nbins = 30)
    edges = calc.edges[1]
    s = step(edges)
    theory = [integrate(x -> chord_length_theory(x, R, λ), n:0.05:n+s)
              for n in 0:s:s*(length(edges) - 2)]

    # Compare cummulative distributions instead of probability
    # densities because density is almost zero for big radii.
    calc_cdf = scan(calc.weights)
    theory_cdf = scan(theory)

    err = relerr.(calc_cdf, theory_cdf)
    @test maximum(err) < 0.25
    @test relerr(mc, mean_chord_length(R, λ)) < 0.25
end
