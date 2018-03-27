using Distributions

"""
invariant_embedding(;dist = Distributions.Normal(), npts::Int = 100, E::Int = 3, tau::Int = 1)

Draw a new random time series from `dist` until the embedding formed with the given
embedding dimension `E` and embedding lag `tau` forms an invariant set, i.e. that
the last point of the embedding falls within the convex hull of the previous points.
"""
function invariant_gaussian_embedding(;dist = Distributions.Normal(), covariance = 0.2, npts::Int = 100, E::Int = 3, tau::Int = 1)
   count = 0
   is_invariant = false
   while !is_invariant && count < 100
      count += 1
      # Create source and target time series from the given distribution
      source = rand(dist, npts, 1)
      dest = covariance .* source[1:end] .+ (1.0 - covariance) .* rand(dist, npts, 1)

      # Embed the time series for TE
      embedding = hcat(source[1:end-tau], dest[1:end-tau], dest[1+tau:end])

      # Check if it's invariant. If not, repeat until it ts.
      if SimplexSplitting.invariantset(embedding)
         is_invariant = true
         return embedding
      end
   end
end

"""
   invariantize_embedding(
      embedding::Array{Float64, 2};
      max_point_remove::Int = ceil(Int, size(embedding, 1)*0.1)
      )

If `remove_points = true`, itteratively remove the last point of an embedding until it is invariant. If it is not
possible to render the embedding invariant, return an empty array. The default is to try to
remove a maximum of ~10% of the points of the original embedding before giving up. The
number of points we're allowed to remove can be set by providing the named argument
`max_point_remove`.

If `remove_points = false`, incrementally move last point of the embedding towards the origin until it
"""
function invariantize_embedding(
   embedding::Array{Float64, 2};
   remove_points = false,
   max_point_remove::Int = ceil(Int, size(embedding, 1)*0.05)
   )

   if size(unique(embedding, 1)) < size(embedding)
      warn("Embedding points are not unique. Returning nothing.")
      return nothing
   end
   pts_removed = 0
   is_invariant = false
   percent_moved = 0.0

   embedding_center = mean(embedding, 1)
   lastpoint = embedding[end, :] #(keep track, in case we need to move it)
   direction = embedding_center.' - lastpoint

   while !is_invariant &&
         pts_removed <= max_point_remove
      pts_removed += 1


      #@show size(e, 1), size(unique(e, 1), 1)
      #println()
      if invariantset(embedding[1:(size(embedding, 1) - pts_removed), :])
         is_invariant = true
      else
         percent_moved += 1
         warn("Moved point $percent_moved % towards origin to fit inside convex hull of previous points")
         embedding[end, :] = lastpoint + direction*(percent_moved / 100)
      end
   end

   if is_invariant
      return embedding[1:(size(embedding, 1) - pts_removed), :]
   else
      warn("Could not make embedding invariant. Returning nothing.")
      return nothing
   end
end
