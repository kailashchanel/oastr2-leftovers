#=
Using the DBSCAN algorithm, group the galaxy data into diffeent clusters
=#

import Pkg
Pkg.add(["Clustering", "Cosmology", "DataFrames", "CSV", "Unitful", "PlotlyJS", "JLD", "Distances", "Plots"])

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

"Logs message with timestamp to text file"
function time_log(msg)
    open(string(pwd(), "/output/time_log.txt"), "a") do io
        write(io, string(now(), ": $msg\n"))
    end
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

function find_optimal_clusters(data, name)
    R = dbscan(load_dist_matrix(data, string(pwd(), "/data/$name-distance-matrix.jld")), 2)
    # R = dbscan(retrieve_data_matrix(data), 2)

    println(length(R))
    println(length(getproperty.(R, :counts))) # does not work

    # for i in 1:20
    #     println(R[i])
    # end
    # println(length(counts(R)))
    # println("# of Clusters: $(length(R.clusters))")
    # println("Seeds: $(R.seeds)")
    # println("Counts: $(R.counts)")

    # data.assignments = R[1].assignments
    # generate_plot(data, string(pwd(), "/output/$name-dbscan.html"), length(R[1].assignments), name)
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
    
    find_optimal_clusters(datasets[1], names[1])

    println(string("Logging end time: ", now()))
    time_log("Program ends")
    println("\nProgram complete.")
end


main()