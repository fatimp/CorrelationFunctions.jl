const cfs  = Dict(
    1 => [:s2, :c2],
    2 => [:s2, :c2, :surfsurf, :surfvoid],
    3 => [:s2, :c2, :surfsurf, :surfvoid]
)

const dirs = Dict(
    1 => [:x],
    2 => [:x, :y, :xy, :yx],
    3 => [:x, :y, :z]
)

function test_cf_on_img(img, cf_map, cf_dir)
    directions = dirs[ndims(img)]

    for periodic in [true, false]
        cmap = cf_map(img, true; periodic)
        cdir = cf_dir(img, true; len=size(img, 1), periodic, directions)

        for (direction, data) in cdir
            @test Map.dir_from_map(cmap, direction; periodic) ≈ data
        end
    end
end

function test_all()
    for dim in 1:3
        img = rand(Bool, (10 for _ in 1:dim)...)
        
        for cf in cfs[dim]
            @testset "random $(dim)D image, $(cf)" begin
                test_cf_on_img(img, Map.eval(cf), Directional.eval(cf))
            end
        end
    end
end

test_all()
