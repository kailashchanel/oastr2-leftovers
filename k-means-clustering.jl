# using Pkg
# Pkg.add(["Clustering", "Cosmology", "DataFrames", "CSV", "Unitful", "PlotlyJS", "JLD", "Distances", "Plots"])

using Clustering, Cosmology, DataFrames, Distances, CSV, JLD, Statistics, Unitful, PlotlyJS
import Plots

# Add dist values to the dataset
function add_dist!(data)
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

    data.dist = comoving_radial_dist.(Ref(c), data.Z_HELIO)
end

function find_clusters(data, name)
    gen_plot = true
    gen_silhouette_plot = true
    data_file_path = string("./data/", name, "-", "distance-matrix.jld")
    
    raw_dist_arr = Vector{Float64}()
    println("Stripping units from 'dist' column to create 'raw_dist' column...")
    for n in eachindex(data[!, "dist"])
        append!(raw_dist_arr, Unitful.ustrip(data[!, "dist"][n]))
    end
    data.raw_dist = raw_dist_arr
    
    X = Vector{Float64}()
    Y = Vector{Float64}()
    Z = Vector{Float64}()
    println("Converting spherical coordinates to rectangular coordinates...")
    for n in eachindex(raw_dist_arr)
        val_x = raw_dist_arr[n] * sin(deg2rad(180 - data[!, "DEC"][n])) * cos(deg2rad(data[!, "RA"][n]))
        val_y = raw_dist_arr[n] * sin(deg2rad(180 - data[!, "DEC"][n])) * sin(deg2rad(data[!, "RA"][n]))
        val_z = raw_dist_arr[n] * cos(deg2rad(180 - data[!, "DEC"][n]))

        append!(X, val_x)
        append!(Y, val_y)
        append!(Z, val_z)
    end
    data.X = X
    data.Y = Y
    data.Z = Z

    println("Converting Vectors to 3-row coordinate Matrix...")
    data_matrix = Matrix{Float64}(undef, 3, 0)
    for n in eachindex(data[!, "RA"])
        A = [data[!, "X"][n]; data[!, "Y"][n]; data[!, "Z"][n];]
        data_matrix = hcat(data_matrix, A)
    end
    
    println(size(data_matrix))
    
    if !(isfile(data_file_path))
        println("Creating distance matrix...")
        P = pairwise(Euclidean(), data_matrix, dims=2)
        println("Saving distance matrix to '", data_file_path, "' for future use...")
        save(data_file_path, "data", P)
    else
        println(string("Loading distance matrix from '", data_file_path, "'..."))
        P = load(data_file_path, "data")
    end
    
    avg_silhouette_scores = Dict()
    for i in 50:50
        R = kmeans(data_matrix, i; maxiter=200, display=:iter)
        # println(typeof(R))
        println(string("# of clusters: ", i))
        println(string("# in each cluster: ", counts(R)))
        # data.assignments = assignments(R)
        
        println("Calculating silhouette scores...")
        S = silhouettes(R, P)

        subclusters = Dict{Integer, Vector}()
        for n in 1:i
            subclusters[n] = Vector{Float64}()
        end

        for n in eachindex(S)
            append!(subclusters[assignments(R)[n]], S[n])
        end

        sub_averages = Dict{Integer, Float64}()
        for (group_num, value) in subclusters
            sub_averages[group_num] = mean(value)
        end

        avg_silhouette_scores[string(i)] = mean(S)
        avg_silhouette_scores[string(i, "-subcluster-averages")] = sub_averages
        avg_silhouette_scores[string(i, "-subcluster-counts")] = counts(R)

        y = Vector{Float64}()
        for (key, val) in sub_averages
            append!(y, val)
        end
        
        if (gen_silhouette_plot)
            println("Generating avg silhouette per subcluster bar chart...")
            p = Plots.bar(avg_silhouette_scores[string(i, "-subcluster-counts")], y, xlabel="# in subcluster", ylabel="avg silhouette score")
            Plots.savefig(p, string("./output/", name, "-", i, "-subcluster-averages.png"))
        end

        println(string("Average Silhouette Score: ", avg_silhouette_scores[string(i)]))

        data.assignments = assignments(R)
        output_file_path = string("./output/", name, "-", i, "clusters-3dspace.html")
        if (!(isfile(output_file_path)) && gen_plot)
            println("Generating 3D interactive scatter plot...")
            generate_plot(data, output_file_path, i, name)
        end
    end
    println("Average Silhouette Scores: ")
    for (key, val) in avg_silhouette_scores
        println(string(key, " => ", val))
    end
end

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

function main()
    #define cosmological model. For this example I will use the Planck 2015 
    #cosmological parameters but this can be easily modified. 

    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

    age(c, 0)
    age(c, 1.42) #z = 1.42

    data1 = CSV.read("GAMA_CZ5Unj.csv", DataFrame)

    # data1.age = age.(Ref(c), data1.Z_HELIO)

    println("Calculating radial distance...")
    add_dist!(data1)

    println("Calculating data groups...")
    G02 = data1[((data1[!, "RA"].<38.8) .& (data1[!, "RA"].>30.2) .& (data1[!, "DEC"].<-3.72) .& (data1[!, "DEC"].>-10.25)),:]
    G09 = data1[((data1[!, "RA"].<141.0) .& (data1[!, "RA"].>129.0) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
    G12 = data1[((data1[!, "RA"].<186.0) .& (data1[!, "RA"].>174.0) .& (data1[!, "DEC"].<2.0) .& (data1[!, "DEC"].>-3.0)),:]
    G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
    G23 = data1[((data1[!, "RA"].<351.9) .& (data1[!, "RA"].>338.1) .& (data1[!, "DEC"].<-30.0) .& (data1[!, "DEC"].> -35.0)),:]

    println("Generating graphs...")
    find_clusters(G02, "G02")
    # generate_graph(G09, "g09-radial.html", "G09")
    # generate_graph(G12, "g12-radial.html", "G12")
    # generate_graph(G15, "g15-radial.html", "G15")
    # generate_graph(G23, "g23-radial.html", "G23")

    println("Program complete.")
end


main()