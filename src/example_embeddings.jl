using Distributions, SimplexSplitting

"""
invariant_embedding(;dist = Distributions.Normal(), npts::Int = 100, E::Int = 3, tau::Int = 1)

Draw a new random time series from `dist` until the embedding formed with the given
embedding dimension `E` and embedding lag `tau` forms an invariant set, i.e. that
the last point of the embedding falls within the convex hull of the previous points.
"""
function inv_gauss(;dist = Distributions.Normal(), cov = 0.2, npts::Int = 100, E::Int = 3, tau::Int = 1)
   count = 0
   is_invariant = false

   while !is_invariant && count < 100
      count += 1
      # Create source and target time series from the given distribution
      source = rand(dist, npts, 1)
      dest = cov .* source[1:end] .+ (1.0 - cov) .* rand(dist, npts, 1)

      # Embed the time series for TE
      embedding = hcat(source[1:end-tau], dest[1:end-tau], dest[1+tau:end])

      # Check if it's invariant. If not, repeat until it ts.
      if SimplexSplitting.invariantset(embedding)
         is_invariant = true
         return embedding
      end
   end
end

include("invariantize_embedding.jl")
"""
invariant_embedding(;dist = Distributions.Normal(), npts::Int = 100, E::Int = 3, tau::Int = 1)

Draw a new random time series from `dist` until the embedding formed with the given
embedding dimension `E` and embedding lag `tau` forms an invariant set, i.e. that
the last point of the embedding falls within the convex hull of the previous points.
"""
function invariant_gaussian_embedding(;npts=50, cov=0.3, E=3, tau=1)
   ig = nothing
   while ig == nothing
      ig = inv_gauss(npts=npts, cov=cov, E=E, tau=tau)
      ig = invariantize_embedding(ig)
   end

   return ig
end
