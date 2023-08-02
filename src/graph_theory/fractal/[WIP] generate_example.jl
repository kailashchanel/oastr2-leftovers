function generate_sierpinski_triangle_adj_matrix(levels)
    n = 2^(levels - 1) * 3
    adj_matrix = zeros(Int, n, n)

    function fill_sierpinski_triangle(adj_matrix, level, top_row, left_col, right_col)
        if level == 1
            adj_matrix[top_row, left_col] = 1
            adj_matrix[top_row, right_col] = 1
            adj_matrix[left_col, right_col] = 1
        else
            mid_col = (left_col + right_col) ÷ 2
            fill_sierpinski_triangle(adj_matrix, level - 1, top_row, left_col, mid_col)
            fill_sierpinski_triangle(adj_matrix, level - 1, top_row, mid_col, right_col)
            fill_sierpinski_triangle(adj_matrix, level - 1, top_row, left_col, right_col)
        end
    end

    fill_sierpinski_triangle(adj_matrix, levels, 1, 2, n)

    return adj_matrix
end

levels = 5  
adj_matrix = generate_sierpinski_triangle_adj_matrix(levels)

println("Sierpinski Triangle Adjacency Matrix:")
println(adj_matrix)