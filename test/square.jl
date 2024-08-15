function make_square(side, ϕ)
    array = centered(falses(side, side))
    for idx in CartesianIndices(array)
	x, y = idx[1]/side, idx[2]/side
	x1 = x*cos(ϕ) + y*sin(ϕ)
	y1 = y*cos(ϕ) - x*sin(ϕ)
	if abs(x1) < 0.2 && abs(y1) < 0.2
	    array[idx] = true
	end
    end

    return array.parent
end

#=
Rationale behind this test:

Let I be an binary image (array) with central symmetry and F_{ss}(I)
be a surface-surface correlation map of I. Let T be a unitary
operator. Then F_{ss}(T I) = T[F_{ss}(I)].

In this test we rotate a square and check a part of an area where
F_{ss} is defined and non-zero. No matter how the square is rotated,
the value of F_{ss} in this area must be equal to 2.
=#
@testset "Check rotation stability of edge detection filters" begin
    for ϕ in 0:0.05:(pi/4)
        square = make_square(3000, ϕ)
        ss = M.surf2(square, true; mode = U.Periodic()) * 3000^2
        slice = ss[60:80, 120:140]
        # F_{ss} for a square is always 2 where defined and non-zero
        @test U.lowfreq_energy_ratio(square) > 0.97
        @test maximum(relerr.(slice, 2)) < 0.06
    end
end
