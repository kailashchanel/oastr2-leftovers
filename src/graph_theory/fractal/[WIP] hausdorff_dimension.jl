using LinearAlgebra

function hausdorff_dimension(adj_matrix)
    eigenvalues = eigvals(adj_matrix)
    max_eigenvalue = maximum(eigenvalues)
    hausdorff_dim = log(count(e -> e >= max_eigenvalue / 2, eigenvalues)) / log(2.0)
    return hausdorff_dim
end

file = CSV.read("adjacency_matrix.csv", DataFrame)
adj_matrix = Matrix(file)
dimension = hausdorff_dimension(adj_matrix)

println("Hausdorff Dimension: ", dimension)