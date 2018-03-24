using InvariantDistribution

@testset  "Markov matrix: discrete, dense" begin
    @test all(sum(test_mm_discrete_bullmap_cube(n_randpts = 10), 2) .≈ 1)
end


#
# using BenchmarkTools
# all(sum(test_mm_discrete_bullmap_cube(n_randpts = 5, prefilter = true, use_orientations = false, split_factor = 2), 2) .≈ 1)
# all(sum(test_mm_discrete_bullmap_cube(n_randpts = 5, prefilter = true, use_orientations = true, split_factor = 2), 2) .≈ 1)
# all(sum(test_mm_discrete_bullmap_cube(n_randpts = 5, prefilter = false, use_orientations = false, split_factor = 2), 2) .≈ 1)
# all(sum(test_mm_discrete_bullmap_cube(n_randpts = 5, prefilter = false, use_orientations = true, split_factor = 2), 2) .≈ 1)
#
# split = 6
# n_pts = 1000
# t1 = @btime all(sum(test_mm_discrete_bullmap_cube(n_randpts = n_pts, prefilter = true, use_orientations = false, split_factor = split), 2) .≈ 1)
# t2 = @btime all(sum(test_mm_discrete_bullmap_cube(n_randpts = n_pts, prefilter = true, use_orientations = true, split_factor = split), 2) .≈ 1)
# #t3 = @btime all(sum(test_mm_discrete_bullmap_cube(n_randpts = n_pts, prefilter = false, use_orientations = false, split_factor = split), 2) .≈ 1)
# #t4 = @btime all(sum(test_mm_discrete_bullmap_cube(n_randpts = n_pts, prefilter = false, use_orientations = true, split_factor = split), 2) .≈ 1)
#
# size(invariant_bullmap_on_cube(split_factor = split).simplex_inds, 1)
# @time all(sum(test_mm_discrete_bullmap_cube(n_randpts = n_pts, prefilter = true, use_orientations = false, split_factor = 7), 2) .≈ 1)
