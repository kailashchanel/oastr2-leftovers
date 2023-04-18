#=
A file to compare MSTs. Used for
comparing the non-brute force algorithm vs the brute force algorithm
for calculating distance matrix.
=#

import Pkg
Pkg.add(["Graphs", "SimpleWeightedGraphs"])
using Graphs, SimpleWeightedGraphs


function equal_edge_weights(graph1, graph2)
    graph1_weights = []
    graph2_weights = []

    for i in range(2, nv(graph1))
        push!(graph1_weights, get_weight(graph1, i - 1, i))
        push!(graph2_weights, get_weight(graph2, i - 1, i))
    end

    sort!(graph1_weights)
    sort!(graph2_weights)

    for i in eachindex(graph1_weights)
        if graph1_weights[i] - graph2_weights[i] > 0.000001  # set this to a very low number because floating point accuracy
            return false;
        end
    end

    return true;
end


function main()
    graph_filename_1 = ""  # FILL THIS IN
    graph_filename_2 = ""  # FILL THIS IN

    graph1 = loadgraph(graph_filename_1, SWGFormat())
    graph2 = loadgraph(graph_filename_2, SWGFormat())

    if nv(graph1) == nv(graph2) && ne(graph1) == ne(graph2)
        println("Graphs are of equal size")
    else
        println("Graphs are not of equal size. They are not the same graph.")
        return
    end

    if equal_edge_weights(graph1, graph2)
        println("Graph has equal edge weights")
    else
        println("Graphs do not have equal weights. They are not the same graph.")
        return
    end

    println("They are likely (but not guaranteed because of some edge cases) the same graph")
end


main()
