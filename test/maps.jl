const cfs  = Dict(
    1 => [:s2, :c2],
    2 => [:s2, :c2, :surf2, :surfvoid],
    3 => [:s2, :c2, :surf2, :surfvoid]
)

const dirs = Dict(
    1 => [D.DirX()],
    2 => [D.DirX(), D.DirY(), D.DirXY(), D.DirYX()],
    3 => [D.DirX(), D.DirY(), D.DirZ()]
)

function test_cf_on_img(img, cf_map, cf_dir)
    directions = dirs[ndims(img)]
    mask = rand(Bool, size(img))

    for mode in (U.Periodic(), U.NonPeriodic(), U.Mask(mask))
        cf(dir) = dir => cf_dir(img, true, dir; len=size(img, 1), mode)
        cmap = cf_map(img, true; mode)
        cdir = cf.(directions)

        for (direction, data) in cdir
            @test all(isapprox_or_nan.(M.dir_from_map(cmap, direction; mode), data))
        end
    end
end

function test_all()
    for dim in 1:3
        img = rand(Bool, (10 for _ in 1:dim)...)
        
        for cf in cfs[dim]
            @testset "random $(dim)D image, $(cf)" begin
                test_cf_on_img(img, M.eval(cf), D.eval(cf))
            end
        end
    end
end

test_all()
