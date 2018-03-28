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
