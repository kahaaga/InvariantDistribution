using Distributions
using InvariantDistribution
using SimplexSplitting


"""
Approximate the map by a Markov matrix obtained by projecting randomly
"""
function mm_discrete_dense_test(;
        #t::SimplexSplitting.Triangulation;
        n_randpts::Int = 10,
        dist::Distributions.Distribution = Distributions.Uniform(0, 1))

    # Create example triangulation with 100 starting points
    dim = 3
    # Create bogus triangulation
    n_embedding_points = 1000
    n_simplices = 1000

    embedding = rand(3, n_embedding_points)
    points = embedding[:, 1:end-1].'
    impoints = embedding[:, 2:end].'
    simplex_inds = rand(1:n_embedding_points-1, n_simplices, dim+1)


    # Create a set of convex coefficients to be used for all simplices
    convex_coeffs = rand(dist, n_randpts, dim + 1)
    convex_coeffs .= convex_coeffs ./ sum(convex_coeffs, 2)
    projected_point = zeros(Float64, 1, dim) # pre-allocate

    occurences = zeros(Int, n_simplices, n_simplices)

    for j in 1:n_simplices
        image_simplex = view(impoints, view(simplex_inds, j, :), :)
        for n in 1:n_randpts
            projected_point .= convex_coeffs[n, :].' * image_simplex

            for i in 1:n_simplices
                simplex = view(points, view(simplex_inds, i, :), :)
                if contains_point(simplex, projected_point)
                    occurences[j, i] += 1
                    break
                end
            end

        end
    end

    return occurences / n_randpts
end


"""
Approximate the map by a Markov matrix obtained by projecting randomly
"""
function mm_discrete_dense(
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
        image_simplex = view(t.points, view(t.simplex_inds, i, :), :)

        for k in 1:n_randpts
            @show convex_coeffs[k, :]
            projected_point = convex_coeffs[k, :].' * image_simplex

            for j in 1:n_simplices
                simplex = view(t.points, view(t.simplex_inds, j, :), :)

                if contains_point(simplex, projected_point)
                    #println("The ", k, "th random point lies inside simplex #", j)
                    occurences[i, j] += 1
                    break
                end
            end

        end
    end

    return occurences / n_randpts
end

"""
Approximate the map by a Markov matrix obtained by projecting randomly
"""
function mm_discrete_dense_prefilter(
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
Approximate the map by a Markov matrix obtained by projecting randomly
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

    # Compute the first sign
    replaced_simplex = simplex
    replaced_simplex[1, :] = point
    signs[1] = sign(det([replaced_simplex ones(4)]))

    # Check remaining signs and stop if sign changes
    for i = 2:size(simplex, 1)
        replaced_simplex = simplex
        replaced_simplex[i, :] = point
        signs[i] = sign(det([replaced_simplex ones(4)]))
        if !(signs[i] == signs[i-1])
            return false
        end
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
        replaced_simplex = simplex
        replaced_simplex[i, :] = point
        α = det([replaced_simplex ones(4)]) * orientation

        if α < 0; return false; end
    end
    return true
end


function childpoint(parentsimplex::Array{Float64, 2})
    dim = size(parentsimplex, 2)
    # Random linear combination coefficients
    R = rand(1, dim + 1)

    # Normalise the coefficients so that they sum to one. We can then create the new point
    # as a convex linear combination of the vertices of the parent simplex.
    normalised_coeffs = (1 ./ sum(R, 2)) .* R

    normalised_coeffs * parentsimplex
end


function outsidepoint(parentsimplex::Array{Float64, 2})
    dim = size(parentsimplex, 2)
    # Random linear combination coefficients
    R = rand(1, dim + 1)

    # Normalise the coefficients so that they sum to one. We can then create the new point
    # as a convex linear combination of the vertices of the parent simplex.
    normalised_coeffs = (1 ./ sum(R, 2)) .* R
    normalised_coeffs[1] += 1

    normalised_coeffs * parentsimplex
end


function test_mm_discrete_map(t::Triangulation;
        n_randpts = 10,
        prefilter = true,
        use_orientations = true,
        split_factor = 3,
        dim = 3)

    if prefilter
        mm_discrete_dense_prefilter(t, n_randpts = n_randpts)
    elseif prefilter & use_orientations
        mm_discrete_dense_prefilter_useorientations(t, n_randpts = n_randpts)
    else
        mm_discrete_dense(t, n_randpts = n_randpts)
    end
end
