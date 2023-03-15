import Pkg
Pkg.add(["Graphs", "SimpleWeightedGraphs", "DataFrames"])

using Graphs, SimpleWeightedGraphs, CSV, DataFrames


function get_df_from_dir(dirname::String) 
    full_df = DataFrame(distance=[], src=[], dst=[])

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
    g = SimpleWeightedGraph(data[!, "src"], data[!, "dst"], data[!, "distance"])
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
    # data = CSV.read("sophia_edge_weights.csv", DataFrame)
    data = get_df_from_dir("C:\\users\\d8amo\\Desktop\\Programming\\Julia\\Astrophysics\\3D MST\\G15")

    graph = generate_graph(data)
    mst = generate_mst(graph)
    savegraph("G15_bruteforce.lgz", mst)
end

main()
