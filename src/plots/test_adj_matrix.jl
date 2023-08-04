using CSV, DataFrames

function is_square(adj_matrix)
    rows, cols = size(adj_matrix)
    return rows == cols
end

adj_matrix = CSV.read("adjacency_matrix.csv", DataFrame)  
if is_square(adj_matrix)
    println("The adjacency matrix is square.")
else
    println("The adjacency matrix is not square.")
end