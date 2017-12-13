function markovmatrix(ts::Vector{Float64}, m::Int64, tau::Int64; print = false)
  const tolerance::Float64 = 1/10^12
  const delta::Float64 = 1/10^5

  embedding = SimplexSplitting.embed(ts, m, tau)

  if SimplexSplitting.invariantset(embedding)
    if print; println("The embedding forms an invariant set"); end
  else
    if print; println("The embedding does not form an invariant set"); end
    return
  end

  # Embed using all points except the last (to allow projection of all vertices
  # in the triangulation).
  triangulation = SimplexSplitting.triangulate(embedding[1:end-1, :])

  triangulation_points::Array{Float64, 2} = triangulation[1]
  indices_simplices::Array{Int, 2} = triangulation[2]

  # Subdivide simplices by their barycentric center of mass
  unique_vertex_inds = unique(indices_simplices) # Indices of unique points in triangulation

  # Shift all indices +1 (to get images)
  indices_images = indices_simplices + 1

  # Joggle coplanar simplices
  joggle_factor = 1/10^4


  n_validsimplices::Int = size(indices_images, 1)
  volume_tolerance::Float64 = delta/n_validsimplices # Mismatch for the markovity of the Markov matrix is at most delta

  # Fill the transition probability matrix by finding the volume of intersection
  # between each simplex in the triangulation and its image (just another
  # simplex by projecting each vertex of the original simplex one time step
  # ahead.
  P = zeros(n_validsimplices, n_validsimplices)
  if print; println("Attempting to create transition probability matrix.\n"); end

  for i = 1:n_validsimplices
    indices_of_vertices_forming_the_image_simplex = indices_images[i, :]
    image = embedding[indices_of_vertices_forming_the_image_simplex, :].'

    # Volume of the image of the simplex
    image_volume = abs(det([ones(1, m+1); image]))

    # test if
    if image_volume < tolerance
      # Loop over image vertices
      for col in 1:size(image, 2)

        # Vertices in triangulation
        for vertex in unique_vertex_inds
          v = triangulation[1][vertex, :]
          if isapprox(image[:, col], v)
            if print; println("image simplex #", i, " = ", image, " (image_vol ", image_volume, ")\n shares vertex with point #", vertex, " in triangulation (coordinates = ", triangulation[1][vertex, :], ")\n"); end
          end
        end
      end
      if print; println("--------------------\n"); end
    end
    for j = 1:n_validsimplices
      indices_of_vertices_forming_the_simplex = indices_simplices[j, :]
      simplex = embedding[indices_of_vertices_forming_the_simplex, :].'

      simplex_volume = abs(det([ones(1, m +1); simplex]))
      if (simplex_volume * image_volume > 0 &&
         simplex_volume/image_volume > volume_tolerance)
        P[i,j] = SimplexIntersection.simplexintersection(simplex, image)/image_volume
      end
    end
  end
  #return P
  return StateSpace(ts, embedding, m, tau,
                    triangulation[1], triangulation[2],
                    indices_simplices, indices_images, P)
end


function markovmatrix_parallel(ts::Array{Float64, 1}, m::Int64, tau::Int64; print = false)
  tolerance = 1/10^15
  delta = 1/10^15

  embedding = embed(ts, m, tau)

  if invariantset(embedding)
    #println("The embedding forms an invariant set")
  else
    #println("The embedding does not form an invariant set")
    return
  end

  # Embed using all points except the last (to allow projection of all vertices
  # in the triangulation).
  triangulation = triangulate(embedding[1:end-1, :])

  triangulation_points = triangulation[1]
  indices_simplices = triangulation[2]


  # Subdivide simplices by their barycentric center of mass
  unique_vertex_inds = unique(indices_simplices) # Indices of unique points in triangulation

  # Shift all indices +1 (to get images)
  indices_images = indices_simplices + 1

  # Joggle coplanar simplices
  joggle_factor = 1/10^4


  n_validsimplices = size(indices_images, 1)
  volume_tolerance = delta/n_validsimplices # Mismatch for the markovity of the Markov matrix is at most delta

  # Fill the transition probability matrix by finding the volume of intersection
  # between each simplex in the triangulation and its image (just another
  # simplex by projecting each vertex of the original simplex one time step
  # ahead.
  #P = SharedArray{Float64}(n_validsimplices, n_validsimplices)
  #println("Attempting to create transition probability matrix.\n")
  P = SharedArray{Float64}(n_validsimplices, n_validsimplices)

  @parallel for i in 1:n_validsimplices
    indices_of_vertices_forming_the_image_simplex = indices_images[i, :]
    image = embedding[indices_of_vertices_forming_the_image_simplex, :].'

    # Volume of the image of the simplex
    image_volume = abs(det([ones(1, m+1); image]))

    for j in 1:n_validsimplices
      indices_of_vertices_forming_the_simplex = indices_simplices[j, :]
      simplex = embedding[indices_of_vertices_forming_the_simplex, :].'

      simplex_volume = abs(det([ones(1, m +1); simplex]))
      if simplex_volume * image_volume > 0 &&
        (simplex_volume/image_volume) > volume_tolerance
        P[i,j] = SimplexIntersection.simplexintersection(simplex, image)/image_volume
      end
    end
  end

  fetch(P)

  return StateSpace(ts, embedding, m, tau,
                    triangulation[1], triangulation[2],
                    indices_simplices, indices_images, P)
end
