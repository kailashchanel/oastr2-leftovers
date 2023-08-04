
using DataFrames, SparseArrays, CSV, Tables

function dataframe_to_adjacency_matrix(df::DataFrame)
    nodes = unique(union(df.src, df.dst))
    
    node_indices = Dict(node => index for (index, node) in enumerate(nodes))
    println("initializing adjacency matrix")
    
    num_nodes = length(nodes)
    adjacency_matrix = spzeros(num_nodes, num_nodes)
    
    total_edges = nrow(df)

    for (i, row) in enumerate(eachrow(df))
        src_index = node_indices[row.src]
        dst_index = node_indices[row.dst]
        distance = row.distance
        
        adjacency_matrix[src_index, dst_index] = distance
        adjacency_matrix[dst_index, src_index] = distance 

        @show i, total_edges 
    end
    
    println("COMPLETE[1]")
    
    return adjacency_matrix, nodes
end

function save_adjacency_matrix_to_csv(adj_matrix, nodes, filename)
    println("saving adjacency matrix as file")
    df = DataFrame(adj_matrix, :auto)
    rename!(df, Symbol.(nodes))
    CSV.write(filename, df)
    
    println("COMPLETE[2]")
end

function main()
    println("Loading data")
    data = CSV.read("/Volumes/FMRIDATA/Sophie/G09_bruteforce.lgz", DataFrame)
    df = data[:, 1:3]
    df_new = rename(df, :1 => "src", :2 => "dst", :3 => "distance")

    dataframe_to_adjacency_matrix(df_new)

    adj_matrix, nodes = dataframe_to_adjacency_matrix(df_new)
    
    filename = "adjacency_matrix.csv"
    save_adjacency_matrix_to_csv(Matrix(adj_matrix), nodes, filename)

    println("Adjacency Matrix has been saved to $filename.")
end

main()