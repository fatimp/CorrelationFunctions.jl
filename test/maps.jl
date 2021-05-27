cfs_map = [:s2, :l2, :c2, :ss, :sv]
cfs_dir = [:s2, :l2, :c2, :surfsurf, :surfvoid]
dirs = Dict(
    1 => [:x],
    2 => [:x, :y, :xy_main, :xy_anti],
    3 => [:x, :y, :z]
)


function test_cf_on_img(img, cf_map, cf_dir)
    for periodic in [true, false]
        cmap = Map.corr_function_map(img, cf_map; periodic)
        cdir = cf_dir(img, true; len=size(img, 1), periodic, directions=dirs[ndims(img)])

        for (direction, data) in cdir
            @test Map.dir_from_map(cmap, direction) ≈ data
        end
    end
end


function test_all()
    dims = [1, 2, 3]
    n = 10

    for dim in dims
        ns = Tuple(ones(Int, dim) * n)
        img = rand(Bool, ns)
        
        for (cf_map, cf_dir) in zip(cfs_map, cfs_dir)
            @testset "random $(dim)D image, $(cf_map)" begin
                test_cf_on_img(img, Map.eval(cf_map), Directional.eval(cf_dir))
            end
        end
    end
end


test_all()
