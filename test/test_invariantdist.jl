
function test_InvariantDistribution(; parallel = false, print = false)
  ts3::Vector{Float64} = [1, 5, 4, 2, 3.2, -5, -6, -7, -4, 4, 5, -1, -2, 3, 2, 1, 5, 4, 2, 3.2,
  -5, -6, -7, -4, 4, 5, -1, -2, 3, 2, 1, 5, 4, 2, 3.2, -5, -6, -7, -4, 4, 5,
  -1, -2, 3, 2, 1, 5, 4, 2, 3.2, -5, -6, -7, -4, 4, 5, -1, -2, 3, 2]

  if parallel
    SS = markovmatrix_parallel(ts3, 3, 2, print = print)
  elseif !parallel
    SS = markovmatrix(ts3, 3, 2, print = print)
  end

  invmeasure = invariantdist(SS.P)

  # Invariant distribution
  invdist = invmeasure[1]

  # Reference to the simplices with stritly positive measure (column vector of
  # indices referring to the rows of the triangulation, each corresponding to
  # a simplex).
  positive_measure_simplex_inds = invmeasure[2]

  # Get the simplices with strictly positive measure
  simplex_inds = SS.triangulation_indices[positive_measure_simplex_inds, :]

  return(invdist)
end

@test sum(test_InvariantDistribution()) ≈ 1.0 atol = 1/10^7
@time test_InvariantDistribution()
