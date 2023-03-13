import Pkg
Pkg.add(["Graphs", "SimpleWeightedGraphs", "DataFrames"])

using Graphs, SimpleWeightedGraphs, CSV, DataFrames


function generate_graph()
    println("Loading data")
    data = CSV.read("sophia_edge_weights.csv", DataFrame)
    println("Data loaded")

    println("Generating graph")
    g = SimpleWeightedGraph(data[!, "src"], data[!, "dst"], data[!, "distance"])
    println("Graph generated")
    
    return g
end


function generate_mst(graph)
    println("Generating MST")
    mst = kruskal_mst(graph)
    println("MST generated")

    println(length(mst))

    println("Generating graph dataframe")
    graph_df = DataFrame(src=[], dst=[], weight=[])

    for item in mst
        push!(graph_df, [item.src, item.dst, item.weight])
    end

    println("Graph dataframe generated")

    println("Generating MST graph")
    mst_graph = SimpleWeightedGraph(graph_df[!, "src"], graph_df[!, "dst"], graph_df[!, "weight"])
    println("Generated MST graph")
    
    return mst_graph
end

function main()
    graph = generate_graph()
    mst = generate_mst(graph)
    savegraph("G15.lgz", mst)
end

main()