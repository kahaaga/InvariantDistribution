
type StateSpace
  timeseries::Array{Float64, 1}
  embedding::Array{Float64, 2}
  embeddingdim::Int64
  embeddinglag::Int64
  triangulation_vertices::Array{Float64, 2}
  triangulation_indices::Array{Int64, 2}
  indices_simplices::Array{Int64} # Indices of simplices forming the triangulation of the state space
  indices_images::Array{Int64}   # Indices of the images of the simplices under the map
  P::Array{Float64, 2}   # Transition probability matrix
end


"""
 Extract a simplex from a triangulated state space instance at the index specified by 'index'
"""
function getsimplex(statespace, index)
  statespace.embedding[statespace.indices_simplices[index, :], :].'
end
