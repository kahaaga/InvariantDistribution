"""
Compute the invariant distribution from a Markov matrix.
"""
function invariantdist(markovmatrix; N = 100, tolerance = 1/10^5, delta = 1/10^5)
    count = 1

    M = size(markovmatrix, 1)
    # Distribution we start with (a random distribution). Normalise it.
    zerodistribution = rand(1, M)
    zerodistribution = zerodistribution ./ sum(zerodistribution, 2)

    distribution = zerodistribution * markovmatrix

    distance = norm(distribution - zerodistribution) / norm(zerodistribution)

    check = floor(Int, 1/delta)

    check_points = floor.(Int, collect(1:N).' ./ check) .* collect(1:N).'

    check_points = check_points[check_points .> 0]
    num_checkpoints = size(check_points, 1)

    check_points_counter = 1

    while count <= N && distance >= tolerance
        count = count + 1
        zerodistribution = distribution

        # Apply the Markov matrix to the current state of the distribution
        distribution = zerodistribution * markovmatrix

        if check_points_counter <= num_checkpoints && count == check_points[check_points_counter]
            check_points_counter = check_points_counter + 1
            colsum_distribution = sum(distribution, 2)[1]
            if abs(colsum_distribution - 1) > delta
                distribution = distribution ./ colsum_distribution
            end
        end

        distance = norm(distribution - zerodistribution) / norm(zerodistribution)
    end

    # Do the last normalisation and check
    colsum_distribution = sum(distribution, 2)[1]

    if abs(colsum_distribution - 1) > delta
        distribution = distribution ./ colsum_distribution
    end

    # Find simplices with strictly positive measure
    finaltriang_simplex_indices = round.(Int, SimplexSplitting.heaviside(distribution) .* collect(1:M).')
    finaltriang_simplex_indices = finaltriang_simplex_indices[finaltriang_simplex_indices .> 0]

    #return vec(distribution.')
    # Extract the elements of the invariant measure corresponding to these indices
    invariant_distribution = distribution[finaltriang_simplex_indices]

    return invariant_distribution, finaltriang_simplex_indices
end
