#=
Finds the 10 most optimal number of clusters for the dataset for a human to review.
This program utilizes multithreading, so when the program is calculating the most optimal
clusters, the console output will be unintelligible.
=#

# import Pkg
# Pkg.add(["Clustering", "Cosmology", "DataFrames", "CSV", "Unitful", "PlotlyJS", "JLD", "Distances", "Plots"])

using Clustering, Cosmology, DataFrames, Distances, CSV, JLD, Statistics, Unitful, PlotlyJS, Dates
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
        val_x = data.dist[n] * sin(deg2rad(90 - data[!, "DEC"][n])) * cos(deg2rad(data[!, "RA"][n]))
        val_y = data.dist[n] * sin(deg2rad(90 - data[!, "DEC"][n])) * sin(deg2rad(data[!, "RA"][n]))
        val_z = data.dist[n] * cos(deg2rad(90 - data[!, "DEC"][n]))

        append!(X, val_x)
        append!(Y, val_y)
        append!(Z, val_z)
    end

    data.X = X
    data.Y = Y
    data.Z = Z
end

"Returns 3-row coordinate Matrix from Vectors"
function retrieve_data_matrix(data)
    println("Converting Vectors to 3-row coordinate Matrix...")
    data_matrix = Matrix{Float64}(undef, 3, 0)
    for n in eachindex(data[!, "RA"])
        A = [data[!, "X"][n]; data[!, "Y"][n]; data[!, "Z"][n];]
        data_matrix = hcat(data_matrix, A)
    end

    return data_matrix
end

"Loads the distance matrix if the filepath exists. If not, it creates a distance matrix and writes to the filepath."
function load_dist_matrix(data, filepath::String)
    if !(isfile(filepath))
        data_matrix = retrieve_data_matrix(data)
        println("Creating distance matrix...")
        time_log("Calculating distance matrix")
        P = pairwise(Euclidean(), data_matrix, dims=2)
        println("Saving distance matrix to '", filepath, "' for future use...")
        save(filepath, "data", P)
    else
        println(string("Loading distance matrix from '", filepath, "'..."))
        P = load(filepath, "data")
    end

    return P
end


"Generates a 3D interactive scatter plot using the given data"
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
            scene = attr(
                # aspectmode="cube",
                # aspectratio=attr(x=1, y=1, z=1),
                xaxis_title="x: Distance (Mpc)",
                yaxis_title="y: Distance (Mpc)",
                zaxis_title="z: Distance (Mpc)",
                xaxis=attr(
                    range=[0,7000]
                ),
                yaxis=attr(
                    range=[0, 7000]
                ),
                zaxis=attr(
                    range=[-7000, 0]
                )
            ),
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

"Generate 3d scatter plot for the subcluster with the highest average global silhouette score"
function prepare_3dscatter(data, name, silhouette_means, assignments)
    for (i, (n_clusters, val)) in enumerate(sort(silhouette_means; byvalue=true, rev=true))
        if i == 1
            data.assignments = assignments[n_clusters]
            generate_plot(data, string(pwd(), "/output/$name-$n_clusters.html"), n_clusters, name)

            break
        end
    end
end

function generate_silhouette_plot(sub_counts, sub_averages, n_clusters)
    y = Vector{Float64}()
    for (key, val) in sub_averages
        append!(y, val)
    end
    p = Plots.bar(sub_counts, y, xlabel="# in subcluster", ylabel="avg silhouette score")
    Plots.savefig(p, string(pwd(), "/output/", n_clusters, "-subcluster-averages.png"))
end


"Returns the silhouette score of n clusters in the data"
function calc_silhouette_scores(P, n_clusters::Int)::NamedTuple
    R = kmedoids(P, n_clusters; maxiter=200, display=:iter)
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
    
    return (total_average = mean(S), sub_averages = sub_averages, counts = counts(R), assignments=assignments(R))
end

"Generate list of top 10 global silhouette scores"
function generate_top10_global_scores(silhouette_means, silhouettes, name, min_clusters, max_clusters)
    str_output = "\n"

    for (i, (n_clusters, val)) in enumerate(sort(silhouette_means; byvalue=true, rev=true))
        if i > 10
            break
        end
        
        silhouette_data = silhouettes[n_clusters]
        str_output *= string(n_clusters, " (", string(silhouette_means[n_clusters]), ")") * ": " * string(silhouette_data.sub_averages) * "\n"
    end

    open(string(pwd(), "/output/$name-silhouette-top-avg-scores-$min_clusters-$max_clusters.txt"), "a") do io
        write(io, str_output)
    end
end

"Generates the silhouette score histogram"
function generate_silhouette_score_histogram(name, silhouette_means, min_clusters, max_clusters)
    sorted_silhouette_means = sort(collect(silhouette_means), by = x->x[1])      # sorts Dict by ascending key value
    
    avg_scores = Float64[]
    for (i, val) in sorted_silhouette_means
        append!(avg_scores, val)
    end

    p = Plots.plot(min_clusters:max_clusters, avg_scores, xlabel="Number of clusters (k)", ylabel="Silhouette score", plot_title="$name Average Global Silhouette Scores", legend=false)
    Plots.savefig(p, string(pwd(), "/output/$name-global-silhouette-averages-$min_clusters-$max_clusters.png"))
end

"Logs message with timestamp to text file"
function time_log(msg)
    open(string(pwd(), "/output/time_log.txt"), "a") do io
        write(io, string(now(), ": $msg\n"))
    end
end

"Generate the outputs given the silhouette means and the silhouettes"
function generate_output(name, data, silhouettes, silhouette_means, assignments, min_clusters, max_clusters) 
    prepare_3dscatter(data, name, silhouette_means, assignments)
    generate_top10_global_scores(silhouette_means, silhouettes, name, min_clusters, max_clusters)
    generate_silhouette_score_histogram(name, silhouette_means, min_clusters, max_clusters)
end

function find_optimal_clusters(P, min_clusters::Int, max_clusters::Int)   
    silhouettes = Dict()      # Key: silhouette ID, Value: silhouette data
    silhouette_means = Dict()  # Key: silhouette ID, Value: silhouette mean
    assignments = Dict{Integer, Vector}()      # Key: silhouette ID, Value: assignments Vector

    Threads.@threads for i in min_clusters:max_clusters
        time_log("Beginning calculations for $i clusters")
        println("Calculating for $i clusters...")
        silhouette_data = calc_silhouette_scores(P, i)

        silhouettes[i] = silhouette_data
        silhouette_means[i] = silhouette_data.total_average
        assignments[i] = silhouette_data.assignments
    end

    return silhouettes, silhouette_means, assignments
end

function analyze_clusters(name, data, min_clusters=5, max_clusters=100)
    f1 = string(pwd(), "/data/$name-silhouette-data-$min_clusters-$max_clusters.jld")
    f2 = string(pwd(), "/data/$name-silhouette-means-data-$min_clusters-$max_clusters.jld")
    f3 = string(pwd(), "/data/$name-assignments-data-$min_clusters-$max_clusters.jld")

    if !(isfile(f1) && isfile(f2) && isfile(f3))
        silhouettes, silhouette_means, assignments = find_optimal_clusters(
            load_dist_matrix(data, string(pwd(), "/data/$name-distance-matrix.jld")), 
            min_clusters, max_clusters
        )

        println(string("Saving silhouette data to $f1..."))
        save(f1, "data", silhouettes)
        println(string("Saving silhouette means data to $f2..."))
        save(f2, "data", silhouette_means)
        println(string("Saving assignments data to $f3..."))
        save(f3, "data", assignments)
    else
        println("Loading silhouette data from $f1 and $f2...")
        silhouettes = load(f1, "data")
        silhouette_means = load(f2, "data")
        assignments = load(f3, "data")
    end

    println("Generating output...")
    generate_output(name, data, silhouettes, silhouette_means, assignments, min_clusters, max_clusters)
end

function main()
    #define cosmological model. For this example I will use the Planck 2015 
    #cosmological parameters but this can be easily modified. 

    data = CSV.read(string(pwd(), "/data/GAMA_CZ5Unj.csv"), DataFrame)

    println("Calculating radial distance...")
    add_dist!(cosmology(h=0.7, OmegaM=0.3, OmegaR=0), data)
    println("Adding rectangular coordinates...")
    add_xyz!(data)

    println("Calculating data groups...")
    G02 = data[((data[!, "RA"].<38.8) .& (data[!, "RA"].>30.2) .& (data[!, "DEC"].<-3.72) .& (data[!, "DEC"].>-10.25)),:]
    G09 = data[((data[!, "RA"].<141.0) .& (data[!, "RA"].>129.0) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
    G12 = data[((data[!, "RA"].<186.0) .& (data[!, "RA"].>174.0) .& (data[!, "DEC"].<2.0) .& (data[!, "DEC"].>-3.0)),:]
    G15 = data[((data[!, "RA"].<223.5) .& (data[!, "RA"].>211.5) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
    G23 = data[((data[!, "RA"].<351.9) .& (data[!, "RA"].>338.1) .& (data[!, "DEC"].<-30.0) .& (data[!, "DEC"].> -35.0)),:]
    names = ["G02", "G09", "G12", "G15", "G23"]
    datasets = [G02, G09, G12, G15, G23]

    println(string("Logging start time: ", now()))
    time_log("Program begins")
    
    for i in eachindex(datasets)     # iterating through the datasets is too taxing on system memory (ex. G09's data matrix is ~40 GB) - this is a problem that must be fixed (or solved by running on a lab computer or server).
        time_log("Analyzing $(names[i]) data")
        analyze_clusters(names[i], datasets[i])
    end

    # time_log("Analyzing $(names[1]) data")
    # analyze_clusters(names[1], datasets[1], 20, 35)

    println(string("Logging end time: ", now()))
    time_log("Program ends")
    println("\nProgram complete.")
end


main()