"""
   invariant_embedding(;dist = Distributions.Normal(), npts::Int = 100, E::Int = 3, tau::Int = 1)

Draw a new random time series from `dist` until the embedding formed with the given
embedding dimension `E` and embedding lag `tau` forms an invariant set, i.e. that
the last point of the embedding falls within the convex hull of the previous points.
"""
function invariant_embedding(;dist = Distributions.Normal(),
                             npts::Int = 100, E::Int = 3, tau::Int = 1)

   count = 0
   is_invariant = false
   while !is_invariant && count < 100
      count += 1
      ts = rand(dist, npts)
      e = embedding(ts, E, tau)

      if SimplexSplitting.invariantset(e.embedding)
         is_invariant = true
         return e
      end
   end
end
