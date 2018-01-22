using Distributions
using InvariantDistribution
using SimplexSplitting
using TransferEntropy
using PlotlyJS

function te_for_dist(source, target; binsizes = 1:1:100, te_lag = 1)

    embedding = hcat(source[1:end-te_lag], source[te_lag:end], target[te_lag:end])

    if invariantset(embedding)
      #println("The embedding forms an invariant set. Continuing.")
    else
      #println("The embedding does not form an invariant set. Quitting.")
      return zeros(length(binsizes))
    end

    # Embed using all points except the last (to allow projection of all vertices
    # in the triangulation).
    points, simplex_inds = triangulate(embedding[1:end-1, :])
    image_points = embedding[2:end, :]

    P = markovmatrix(points, image_points, simplex_inds)

    invmeasure, inds_nonzero_simplices = invariantdist(P)
    #centroids = rand(dist, npts, 3)
    #invmeasure = abs.(rand(dist, npts))
    #invmeasure[randperm(400)[50]] = 0
    #invmeasure = invmeasure ./ sum(invmeasure) # normalise to true probability dist.

    TE = zeros(length(binsizes))

    count = 0
    for binsize in binsizes
        count +=1
        te = te_from_triangulation(embedding, invmeasure, binsize)
        TE[count] = te
    end
    return TE
end


function te_for_dist(embedding; binsizes = 1:1:100)

    if invariantset(embedding)
      #println("The embedding forms an invariant set. Continuing.")
    else
      #println("The embedding does not form an invariant set. Quitting.")
      return zeros(length(binsizes))
    end

    # Embed using all points except the last (to allow projection of all vertices
    # in the triangulation).
    points, simplex_inds = triangulate(embedding[1:end-1, :])
    image_points = embedding[2:end, :]

    P = markovmatrix(points, image_points, simplex_inds)

    invmeasure, inds_nonzero_simplices = invariantdist(P)
    #centroids = rand(dist, npts, 3)
    #invmeasure = abs.(rand(dist, npts))
    #invmeasure[randperm(400)[50]] = 0
    #invmeasure = invmeasure ./ sum(invmeasure) # normalise to true probability dist.

    TE = zeros(length(binsizes))

    count = 0
    for binsize in binsizes
        count +=1
        te = te_from_triangulation(embedding, invmeasure, binsize)
        TE[count] = te
    end
    return TE
end



function plot_tes(TEs, expected)
    ####################
    # PLOT THE TE CURVES
    ####################
    data = PlotlyJS.GenericTrace[]
    layout = PlotlyJS.Layout(autosize = true,
                              xaxis = PlotlyJS.attr(title = "# bins", ticks = "outside"),
                              yaxis = PlotlyJS.attr(title = "TE (nats)"))
    for i = 1:size(TEs, 2)
        trace = PlotlyJS.scatter(;x = 1:size(TEs, 1),
                                   y = TEs[1:end, i],
                                   mode = "markers+lines",
                                   showlegend = false,
                                   marker = PlotlyJS.attr(size = 3, color = "black"),
                                   line = PlotlyJS.attr(size = 0.8, color = "white"))
        push!(data, trace)
    end
    push!(data, PlotlyJS.scatter(; x = 1:size(TEs, 1),
                                    y = fill(expected, length(1:size(TEs, 1)))))
    PlotlyJS.plot(data, layout)
end

function boxplot_tes(TEs, expected)
    ####################
    # PLOT THE TE CURVES
    ####################
    data = PlotlyJS.GenericTrace[]
    layout = PlotlyJS.Layout(xaxis = PlotlyJS.attr(title = "# bins", ticks = "outside"),
                              yaxis = PlotlyJS.attr(title = "TE (nats)"))
    for i = 1:size(TEs, 2)
        trace = PlotlyJS.box(;y = TEs[1:end, i])
        push!(data, trace)
    end
    push!(data, PlotlyJS.scatter(; x = 1:size(TEs, 1),
                                    y = fill(expected, length(1:size(TEs, 1)))))
    PlotlyJS.plot(data, layout)
end

# COMPUTE TE
#############


function construct_embedding(; a = 1, b = 1, c = 1, epsilon = 0.3, npts = 500, te_lag = 1)
    #exact_TE =  (1/2)*log((1+epsilon^2*a/b)*(1+epsilon^2*c/b)/(1+(a/b+c/b)*epsilon^2))

    dist = Normal()
    X = rand(dist, npts, 1)
    Y = rand(dist, npts, 1)

    x = X[1:end-te_lag]
    y = Y[1:end-te_lag]
    z = epsilon .* x .+ (1 - epsilon) .* rand(dist, npts - 1, 1)

    return hcat(x, y, z)
end


function rep_te(;epsilon = 0.3, reps = 10, binsizes = 1:1:100)
    TEs = zeros(length(binsizes), reps)
    for i = 1:reps
        println("rep#", i)
        TEs[1:end, i] = te_for_dist(construct_embedding(epsilon = epsilon), binsizes = binsizes)
    end
    return TEs
end

# PLOT TE
###########
reps = 50
#epsilon = 0.3 # covariance in this case
a = b = c = 1

#TEs_03 = rep_te(epsilon = 0.3, reps = reps)
TEs_05 = rep_te(epsilon = 0.5, reps = reps)

exact_TE = 0.5*log(1 + (epsilon^2/(1-epsilon)^2)) #exact_TE =  (1/2)*log((1+epsilon^2*a/b)*(1+epsilon^2*c/b)/(1+(a/b+c/b)*epsilon^2))
#plot_tes(TEs_03, 0.5*log(1 + (epsilon^2/(1-epsilon)^2)))
#plot_tes(TEs_05, 0.5*log(1 + (epsilon^2/(1-epsilon)^2)))
plot_tes(median(TEs_05, 2), (1/2)*log(1 + (0.5^2/(1-0.5)^2)))
boxplot_tes(TEs_03, (1/2)*log(1 + (0.3^2/(1-0.3)^2)))

TEs_03 = TEs
