using Pkg
Pkg.add(["Clustering", "Cosmology", "DataFrames", "CSV", "Unitful", "PlotlyJS", "JLD", "Distances", "Plots"])

using Clustering, Cosmology, DataFrames, Distances, CSV, JLD, Statistics, Unitful, PlotlyJS
import Plots

"Add radial distance values to the dataset"
function add_dist!(c, data)
    data.dist = ustrip(comoving_radial_dist.(Ref(c), data.Z_HELIO))
end


"Adds absolute x, y, and z coordinates"
function add_xyz!(data)
    X = Vector{Float64}()
    Y = Vector{Float64}()
    Z = Vector{Float64}()

    println("Converting spherical coordinates to rectangular coordinates...")
    for n in eachindex(data.dist)
        val_x = data.dist[n] * sin(deg2rad(180 - data[!, "DEC"][n])) * cos(deg2rad(data[!, "RA"][n]))
        val_y = data.dist[n] * sin(deg2rad(180 - data[!, "DEC"][n])) * sin(deg2rad(data[!, "RA"][n]))
        val_z = data.dist[n] * cos(deg2rad(180 - data[!, "DEC"][n]))

        append!(X, val_x)
        append!(Y, val_y)
        append!(Z, val_z)
    end

    data.X = X
    data.Y = Y
    data.Z = Z
end


"Loads the distance matrix if the filepath exists. If not, it creates a distance matrix and writes to the filepath."
function load_dist_matrix(filepath::String)
    if !(isfile(filepath))
        println("Creating distance matrix...")
        P = pairwise(Euclidean(), data_matrix, dims=2)
        println("Saving distance matrix to '", filepath, "' for future use...")
        save(filepath, "data", P)
    else
        println(string("Loading distance matrix from '", filepath, "'..."))
        P = load(filepath, "data")
    end

    return P
end


"Generates a plot using the given data"
function generate_plot(data, output_file_path, n_clusters, name)
    p = plot(
        data,
        x=:X,
        y=:Y, 
        z=:Z,
        color=:assignments,
        type="scatter3d", 
        mode="markers",
        Layout(
            title=string("Galaxies Plotted in 3D Space (", n_clusters, " clusters): ", name), 
            # xaxis=attr(title="Radial Distance"),
            # yaxis=attr(title="RA"),
            # zaxis=attr(title="Dec"),
            font=attr(
                family="Courier New",
                size=18,
            )
        )
    )

    open(output_file_path, "w") do io
        PlotlyBase.to_html(io, p.plot)
    end
end


function generate_silhouette_plot(sub_counts, sub_averages, n_clusters)
    y = Vector{Float64}()
    for (key, val) in sub_averages
        append!(y, val)
    end
    p = Plots.bar(sub_counts, y, xlabel="# in subcluster", ylabel="avg silhouette score")
    Plots.savefig(p, string("./output/", n_clusters, "-subcluster-averages.png"))
end


"Returns the silhouette score of n clusters in the data"
function calc_silhouette_scores(P, data_matrix, n_clusters::Int)::NamedTuple
    R = kmeans(data_matrix, n_clusters; maxiter=200, display=:iter)
    println(string("# of clusters: ", n_clusters))
    println(string("# in each cluster: ", counts(R)))
    
    println("Calculating silhouette scores...")
    S = silhouettes(R, P)

    subclusters = Dict{Integer, Vector}()
    for n in 1:n_clusters
        subclusters[n] = Vector{Float64}()
    end

    for n in eachindex(S)
        append!(subclusters[assignments(R)[n]], S[n])
    end

    sub_averages = Dict{Integer, Float64}()
    for (group_num, value) in subclusters
        sub_averages[group_num] = mean(value)
    end

    return (total_average = mean(S), sub_averages = sub_averages, counts = counts(R))
end


function find_optimal_clusters(data, name, min_clusters::Int, max_clusters::Int)
    println("Converting Vectors to 3-row coordinate Matrix...")
    data_matrix = Matrix{Float64}(undef, 3, 0)
    for n in eachindex(data[!, "RA"])
        A = [data[!, "X"][n]; data[!, "Y"][n]; data[!, "Z"][n];]
        data_matrix = hcat(data_matrix, A)
    end
    
    silhouettes = Dict()      # Key: silhouette ID, Value: silhouette data
    silhoutte_means = Dict()  # Key: silhouette ID, Value: silhouette mean

    P = load_dist_matrix(string("./data/", name, "-", "distance-matrix.jld"))

    for i in min_clusters:max_clusters
        println("Calculating cluster $i")
        silhouette_data = calc_silhouette_scores(P, data_matrix, i)

        silhouettes[i] = silhouette_data
        silhoutte_means[i] = silhouette_data.total_average
    end

    # Come up with the output
    # TODO: Refactor this in the future
    str_output = ""

    for (i, (n_clusters, val)) in enumerate(sort(silhoutte_means; byvalue=true))
        if i > 10
            break
        end
        R = kmeans(data_matrix, n_clusters; maxiter=200, display=:iter)
        data.assignments = assignments(R)
        generate_plot(data, "output/" * string(n_clusters) * ".html", n_clusters, "")
        
        silhouette_data = silhouettes[n_clusters]

        generate_silhouette_plot(silhouette_data.counts, silhouette_data.sub_averages, n_clusters)
        str_output *= string(n_clusters) * ": " * string(silhouette_data.sub_averages) * "\n"
    end

    open("output/silhouette_scores.txt", "w") do io
        write(io, str_output)
    end
end


function main()
    #define cosmological model. For this example I will use the Planck 2015 
    #cosmological parameters but this can be easily modified. 

    data = CSV.read("GAMA_CZ5Unj.csv", DataFrame)

    println("Calculating radial distance...")
    add_dist!(cosmology(h=0.7, OmegaM=0.3, OmegaR=0), data)
    println("Adding coordinates...")
    add_xyz!(data)

    println("Calculating data groups...")
    G02 = data[((data[!, "RA"].<38.8) .& (data[!, "RA"].>30.2) .& (data[!, "DEC"].<-3.72) .& (data[!, "DEC"].>-10.25)),:]
    G09 = data[((data[!, "RA"].<141.0) .& (data[!, "RA"].>129.0) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
    G12 = data[((data[!, "RA"].<186.0) .& (data[!, "RA"].>174.0) .& (data[!, "DEC"].<2.0) .& (data[!, "DEC"].>-3.0)),:]
    G15 = data[((data[!, "RA"].<223.5) .& (data[!, "RA"].>211.5) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
    G23 = data[((data[!, "RA"].<351.9) .& (data[!, "RA"].>338.1) .& (data[!, "DEC"].<-30.0) .& (data[!, "DEC"].> -35.0)),:]

    println("Generating graphs...")
    find_optimal_clusters(G02, "G02", 2, 10)

    println("Program complete.")
end


main()
