# Random array with two phases
rand_array = rand(Float32, (100,100,100))
rand_array = map(x -> (x<0.3) ? 0 : 1, rand_array)

# s_2/l_2^i(x) = 0 for all x and i ∉ {0,1}
@test all(x -> x == 0, ms2(rand_array, 50, 2))
@test all(x -> x == 0, ml2(rand_array, 50, 2))

# Probability of randomly choosing 0 or 1
p = count(x -> x == 0, rand_array) / length(rand_array) # P{dot ∈ phase 0}
q = 1 - p                                               # P{dot ∈ phase 1}
prob = [p, q]

# Correlations functions
s = reduce(hcat, ms2(rand_array, 50, i) for i in 0:1)
l = reduce(hcat, ml2(rand_array, 50, i) for i in 0:1)

# Test probabilities of choosing 0 or 1
for i in 1:2
    @test s[1,i] ≈ l[1,i]
    @test abs(s[1,i] - prob[i]) < 1e-3
end

# s_2^i(2) = l_2^i(2)
for i in 1:2
    @test s[2,i] ≈ l[2,i]
end

# Check that l2 is non-increasing           
for i in 1:2
    @test all(x -> x >= 0, map(-, l[:,i], l[2:end,i]))
end

# S_2^i(x) ≈ prob[i]^2 for all x >= 1
for i in 1:2
    all(x -> abs(x-prob[i]^2) < 1e-3, s[2:end,i])
end
