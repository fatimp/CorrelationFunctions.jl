function gencenters(side, λ)
    n = pois_rand(λ * side^2)
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

function s2_theory(r, R, λ)
    tmp = r/(2R)
    tmp2 = (r > 2R) ? 2 : 2/π*(π + tmp*sqrt(1 - tmp^2) - acos(tmp))
    η = λ * π * R^2
    return exp(-η*abs(tmp2))
end

@testset "S2 on random overlapping disks generated by Poisson process (1)" begin
    R = 20; λ = 5e-4
    f = mean ∘ s2
    spheres = gendisks(7000, R, λ)
    s_calc = f(spheres, 2000, 0)
    s_theory = [s2_theory(r, R, λ) for r in 1:2000]

    err = relerr.(s_calc, s_theory)
    @test maximum(err) < 0.07
end

@testset "S2 on random overlapping disks generated by Poisson process (2)" begin
    R = 50; λ = 1e-4
    f = mean ∘ s2
    spheres = gendisks(7000, R, λ)
    s_calc = f(spheres, 2000, 0)
    s_theory = [s2_theory(r, R, λ) for r in 1:2000]

    err = relerr.(s_calc, s_theory)
    @test maximum(err) < 0.07
end
