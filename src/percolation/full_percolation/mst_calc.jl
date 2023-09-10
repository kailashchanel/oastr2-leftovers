using Graphs, SimpleWeightedGraphs, CSV, DataFrames, JLD2


function get_df_from_dir(dirname::String) 
    full_df = DataFrame(src=[], dst=[], distance=[])

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

    sources::Vector{Int64} = [v isa Int64 ? v : parse(Int64, v) for v in data.src]
    destinations::Vector{Int64} = [v isa Int64 ? v : parse(Int64, v) for v in data.dst]
    weights::Vector{Float64} = [v isa Float64 ? v : parse(Float64, v) for v in data.distance]

    g = SimpleWeightedGraph(sources, destinations, weights)

    println("Graph generated")
    
    return g
end
graph = generate_graph(data)

struct GalaxyDistance
    idx1::Int32
    idx2::Int32
    distance::Float32
end

println("Generating MST")
mst = kruskal_mst(graph)
println("MST generated. Length: ", length(mst))

MST_list = GalaxyDistance[]

for i in 1:length(mst)-1
    println("i = $i")
    for j in mst
        push!(MST_list, GalaxyDistance(edge.src, edge.dst, edge.weight))
    end
end
    
save_object("distance_list_mst.jld2", MST_list)

histogram([d.distance for d in MST_list], bins=100, xlabel="Distance galaxies (Mpc)", ylabel="Count",legend = false)
savefig("mst_list_old.png")

