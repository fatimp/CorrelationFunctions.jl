"""
    gen_checkboard(side)

Generate cube array of ones and zeros with checkboard-like
pattern. The size of checkboard cell is 2.
"""
function gen_checkboard(side)
    checkboard = Array{Int, 3}(undef, side, side, side)
    for i in 1:side
        for j in 1:side
            for k in 1:side
                x = ((i-1) & 2 != 0) ? 1 : 0
                y = ((j-1) & 2 != 0) ? 1 : 0
                z = ((k-1) & 2 != 0) ? 1 : 0
                checkboard[k,j,i] = x ⊻ y ⊻ z
            end
        end
    end
    return checkboard
end

cb = gen_checkboard(300)

# S_2 and L_2 for phases 0 and 1
s = reduce(hcat, pms2(cb, 50, i) for i in 0:1)
l = reduce(hcat, pml2(cb, 50, i) for i in 0:1)

# Phases 0 and 1 are distributed evenly
@test s[:,1] ≈ s[:,2]
@test l[:,1] ≈ l[:,2]

# Check L_2
l = l[:,1]
@test l[1] ≈ 1/2
# Would equal to 1/4 on infinite cube
@test abs(l[2] - 1/4) < 0.05
@test all(x -> x == 0, l[3:end])

# Check S_2
s = s[:,1]
@test all(x -> x ≈ 1/4, s[2:2:end])
@test all(x -> x == 0,  s[3:4:end])
@test all(x -> x ≈ 1/2, s[1:4:end])
