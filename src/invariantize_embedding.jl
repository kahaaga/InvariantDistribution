"""
   invariantize_embedding(
      embedding::Array{Float64, 2};
      max_point_remove::Int = ceil(Int, size(embedding, 1)*0.05)
      )

If `remove_points = true`, iteratively remove the last point of an embedding
until it is invariant. If it is not possible to render the embedding invariant,
return an empty array. The default is to try to remove a maximum of ~5% of the
points of the original embedding before giving up. The number of points we're
allowed to remove can be set by providing the named argument `max_point_remove`.

If `remove_points = false`, incrementally move last point of the embedding
towards the origin until it lies within the convex hull of all preceding points.
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

   #=
   # Keep track of the embedding's centerpoint and the original position of the
   # last point in the embedding, so we can move the last point along a line
   # from its orinal position towards the embedding's center, until the point
   # lies inside the convex hull of the preceding points.
   =#
   embedding_center = mean(embedding, 1)
   lastpoint = embedding[end, :]
   direction = embedding_center.' - lastpoint

   pts_removed = 0
   is_invariant = false
   percent_moved = 0.0

   while !is_invariant && pts_removed <= max_point_remove
      pts_removed += 1

      if invariantset(embedding[1:(size(embedding, 1) - pts_removed), :])
         is_invariant = true
      else
         percent_moved += 1
         warn("Moved point $percent_moved % towards origin to fit inside
               convex hull of previous points")
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
