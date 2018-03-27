"""
Compute the invariant probability distribution from a square Markov matrix `M`.
This is done by repeated application of `M` on an initially random distribution
until the distribution converges.
"""
function estimate_invariant_probs(
        M::AbstractArray{Float64, 2};
        N::Int = 100,
        tolerance::Float64 = 1/10^5,
        delta::Float64 = 1/10^5
        )

    counter = 1

    n_simplices = size(M, 1)

    #=
    # Start with a random distribution `Ρ` (big rho). Normalise it so that it
    # sums to 1 and forms a true probability distribution over the simplices.
    =#
    Ρ = rand(Float64, 1, n_simplices)
    Ρ = Ρ ./ sum(Ρ, 2)

    #=
    # Start estimating the invariant distribution. We could either do this by
    # finding the left-eigenvector of M, or by repeated application of M on Ρ
    # until the distribution converges. Here, we use the latter approach,
    # meaning that we iterate until Ρ doesn't change substantially between
    # iterations.
    =#
    distribution = Ρ * M

    distance = norm(distribution - Ρ) / norm(Ρ)

    check = floor(Int, 1 / delta)
    check_pts = floor.(Int, collect(1:N).' ./ check) .* collect(1:N).'
    check_pts = check_pts[check_pts .> 0]
    num_checkpts = size(check_pts, 1)
    check_pts_counter = 1

    while counter <= N && distance >= tolerance
        counter += 1
        Ρ = distribution

        # Apply the Markov matrix to the current state of the distribution
        distribution = Ρ * M

        if (check_pts_counter <= num_checkpts &&
           counter == check_pts[check_pts_counter])

            check_pts_counter += 1
            colsum_distribution = sum(distribution, 2)[1]
            if abs(colsum_distribution - 1) > delta
                distribution = distribution ./ colsum_distribution
            end
        end

        distance = norm(distribution - Ρ) / norm(Ρ)
    end

    # Do the last normalisation and check
    colsum_distribution = sum(distribution, 2)[1]

    if abs(colsum_distribution - 1) > delta
        distribution = distribution ./ colsum_distribution
    end

    # Find simplices with strictly positive measure.
    simplex_inds_nonzero = heaviside(distribution) .* collect(1:M).'
    simplex_inds_nonzero = round(Int, simplex_inds_nonzero)
    simplex_inds_nonzero = simplex_inds_nonzero[simplex_inds_nonzero .> 0]

    # Extract the elements of the invariant measure corresponding to these indices
    return vec(distribution), simplex_inds_nonzero
end
