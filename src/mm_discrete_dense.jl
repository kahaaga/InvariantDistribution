using Distributions
using InvariantDistribution
using SimplexSplitting

"""
Approximate the map by a Markov matrix.

The argument `sample_how` indicates whether the approximation should be random (according
to a Uniform distribution in the space of barycentric coordinates) or "even" in the space
where the simplex lives (points are generated using shape-preserving splitting of simplices).
"""
function mm_discrete_dense(
        t::SimplexSplitting.Triangulation;
        n_randpts::Int = 100,
        dist::Distributions.Distribution = Distributions.Uniform(0, 1),
        sample_randomly::Bool = false,
        prefilter::Bool = true)


    n_simplices = size(t.simplex_inds, 1)
    dim = size(t.points, 2)

    # Create a set of convex coefficients to be used for all simplices
    if sample_randomly
        convex_coeffs = rand(dist, n_randpts, dim + 1)
        convex_coeffs .= convex_coeffs ./ sum(convex_coeffs, 2)
    else
        minimum_split_factor = ceil(Int, n_randpts^(1 / dim))[1]
        convex_coeffs = SimplexSplitting.even_sampling_rules(dim, minimum_split_factor)
    end

    occurences = zeros(Int, n_simplices, n_simplices)
    projected_point = zeros(Float64, 1, dim) # pre-allocate

    for i in 1:n_simplices
        image_simplex = view(t.impoints, view(t.simplex_inds, i, :), :)

        if prefilter
            inds_potential_simplices = potentially_intersecting_simplices(t, i)
        else
            inds_potential_simplices = 1:n_simplices
        end

        for k in 1:n_randpts
            projected_point = convex_coeffs[k, :].' * image_simplex

            for j in inds_potential_simplices
                simplex = view(t.points, view(t.simplex_inds, j, :), :)

                if contains_point(simplex, projected_point)
                    occurences[i, j] += 1
                    break
                end
            end

        end
    end

    return occurences / n_randpts
end



"""
Approximate the transfer operator a Markov matrix using the information we already
have on the orientations of the simplices.
"""
function mm_discrete_dense_prefilter_useorientations(
    t::SimplexSplitting.Triangulation;
    n_randpts::Int = 100,
    dist::Distributions.Distribution = Distributions.Uniform(0, 1))

    n_simplices = size(t.simplex_inds, 1)
    dim = size(t.points, 2)

    # Create a set of convex coefficients to be used for all simplices
    convex_coeffs = rand(dist, n_randpts, dim + 1)
    convex_coeffs .= convex_coeffs ./ sum(convex_coeffs, 2)
    projected_point = zeros(Float64, 1, dim) # pre-allocate

    occurences = zeros(Int, n_simplices, n_simplices)

    for i in 1:n_simplices
        image_simplex = view(t.impoints, view(t.simplex_inds, i, :), :)

        # Considering the distance between centroids and radii of simplices,
        # which simplices can possibly intersect with this image simplex?
        inds_potential_simplices = potentially_intersecting_simplices(t, i)

        for n in 1:n_randpts
            projected_point .= convex_coeffs[n, :].' * image_simplex

            for j in inds_potential_simplices
                simplex = view(t.points, view(t.simplex_inds, j, :), :)
                if contains_point(simplex, projected_point, t.orientations[j])
                    occurences[i, j] += 1
                    break
                end
            end

        end
    end

    return occurences / n_randpts
end

"""
Given a triangulation and an index `image_i` corresponding to an image simplex,
find the indices of original simplices that might have intersection with that image simplex.
"""
function potentially_intersecting_simplices(t::Triangulation, image_i::Int)
    inds_potential_simplices = Int[]

    n_simplices = length(t.radii)

    for i = 1:n_simplices
        dist_difference = ((t.centroids_im[image_i] - t.centroids[i]).' *
                            (t.centroids_im[image_i] - t.centroids[i]) - (t.radii_im[image_i] + t.radii[i])^2)[1]
        if dist_difference < 0
            push!(inds_potential_simplices, i)
        end
    end
    return inds_potential_simplices
end

"""
Check whether a given `point` (1x(dim+1) array) lies inside a `simplex` ((dim+1) x dim array).
"""
function contains_point(simplex::AbstractArray{Float64, 2}, point::AbstractArray{Float64, 2})
    tolerance = 1/10^10

    # Keep track of signs of determinants. If all have the same sign, the point
    # lies inside the simplex. If not, it lies outside.
    signs = zeros(Float64, 4)
    n_pts = size(simplex, 1)
    # Compute the first sign
    signs[1] = sign(det([vcat(point, simplex[2:end, :]) ones(n_pts)]))

    # Check remaining signs and stop if sign changes
    for i = 2:(n_pts - 1)
        signs[i] = sign(det([vcat(simplex[1:(i - 1), :], point, simplex[(i + 1):end, :]) ones(n_pts)]))
        if !(signs[i] == signs[i-1])
            return false
        end
    end
    signs[end] = sign(det([vcat(simplex[1:(end - 1), :], point) ones(n_pts)]))
    if !(signs[end] == signs[end - 1])
        return false
    end
    return true
end


"""
Check whether a given `point` (1x(dim+1) array) lies inside a `simplex` ((dim+1) x dim array)
given the orientation of the simplex.
"""
function contains_point(
        simplex::AbstractArray{Float64, 2},
        point::AbstractArray{Float64, 2},
        orientation::Float64
        )

    for i = 1:size(simplex, 1)
        replaced_simplex = deeopcopy(simplex)
        replaced_simplex[i, :] = point
        α = det([replaced_simplex ones(4)]) * orientation

        if α < 0; return false; end
    end
    return true
end
