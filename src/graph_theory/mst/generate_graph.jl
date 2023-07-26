#=
Generates a MST from the output of the brute-force distance matrix algorithm.
=#

import Pkg
Pkg.add(["Graphs", "SimpleWeightedGraphs", "DataFrames"])

using Graphs, SimpleWeightedGraphs, CSV, DataFrames


function get_df_from_dir(dirname::String) 
    full_df = DataFrame(distance=[], src=[], dst=[])

    for file in readdir(dirname)
        println("Reading file $file")
        df = CSV.read(string(dirname, "/", file), DataFrame)
        append!(full_df, df)
        println("Read file $file")
    end

    return full_df
end


function generate_graph(data::DataFrame)
    println("Generating graph")

    # If you get an error that "src" cannot be parsed, (or something similar) it's likely that your csv file contains
    # data from multiple runs of data_3d_distance, since data_3d_distance.jl appends instead of writes.
    # Either prune the file or delete and re-run data_3d_distance.jl.
    sources::Vector{Int64} = [v isa Int64 ? v : parse(Int64, v) for v in data.src]
    destinations::Vector{Int64} = [v isa Int64 ? v : parse(Int64, v) for v in data.dst]
    weights::Vector{Float64} = [v isa Float64 ? v : parse(Float64, v) for v in data.distance]

    g = SimpleWeightedGraph(sources, destinations, weights)
    println("Graph generated")
    
    return g
end


function generate_mst(graph)
    println("Generating MST")
    mst = kruskal_mst(graph)
    println("MST generated. Length: ", length(mst))

    src::Array{Int64} = []
    dst::Array{Int64} = []
    weight::Array{Float64} = []

    for edge in mst
        push!(src, edge.src)
        push!(dst, edge.dst)
        push!(weight, edge.weight)
    end

    println("Generating MST graph")
    mst_graph = SimpleWeightedGraph(src, dst, weight)
    println("Generated MST graph")
    
    return mst_graph
end

function main()
    # Replace this to the path of a directory that *only* contains the CSV data
    data = get_df_from_dir("C:/Users/d8amo/Desktop/Programming/Julia/Astrophysics/oastr2-leftovers/rand_points")

    graph = generate_graph(data)
    mst = generate_mst(graph)
    savegraph("output.lgz", mst)
end

main()
