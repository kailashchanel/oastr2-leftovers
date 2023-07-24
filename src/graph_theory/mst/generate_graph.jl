#=
Generates a MST from the output of the brute-force distance matrix algorithm.
=#

using Graphs, SimpleWeightedGraphs, CSV, DataFrames


function get_df_from_dir(dirname::String) 
    full_df = DataFrame(dist=[], RA=[], DEC=[])

    for file in readdir(dirname)
        println("Reading file $file")
        df = CSV.read(string(dirname, "\\", file), DataFrame)
        append!(full_df, df)
        println("Read file $file")
    end

    return full_df
end


function generate_graph(data::DataFrame)
    println("Generating graph")
    g = SimpleWeightedGraph(data[!, "RA"], data[!, "DEC"], data[!, "dist"])
    println("Graph generated")
    
    return g
end


function generate_mst(graph)
    println("Generating MST")
    mst = kruskal_mst(graph)
    println("MST generated. Length: ", length(mst))

    RA::Array{Int64} = []
    DEC::Array{Int64} = []
    weight::Array{Float64} = []

    for edge in mst
        push!(RA, edge.RA)
        push!(DEC, edge.DEC)
        push!(weight, edge.weight)
    end

    println("Generating MST graph")
    mst_graph = SimpleWeightedGraph(RA, DEC, weight)
    println("Generated MST graph")
    
    return mst_graph
end

# data = get_df_from_dir("/Volumes/FMRIDATA/output09.csv")
data = CSV.read("/Volumes/FMRIDATA/output09.csv", DataFrame, types = Dict(:src => Int64, :dst => Int64, :distance => Float32))
graph = generate_graph(data)
data = 0
mst = generate_mst(graph)
savegraph("G09_bruteforce.lgz", mst)



