using LightGraphs
using SimpleGraphs
using Random

function percolate_mst(graph::SimpleGraph, p::Float64)
    rng = Random.default_rng()
    
    percolated_mst = SimpleGraph(graph)

    # randomly remove edges with probability p
    for edge in edges(percolated_mst)
        if rng.uniform() < p
            rem_edge!(percolated_mst, edge)
        end
    end

#check if system percolates
    visited = depth_first_search(percolated_mst, 1)
    for i in 2:nv(graph)
        if !(i in visited)
            return false
        end
    end
    return true
end

#adj is the adjacency matrix
graph = SimpleGraph(adj)
mst = kruskal_mst(graph)
p = 0.5
trials = 1000
percolation_probability = percolate_probability(mst, p, trials)

function percolate_probability(graph::SimpleGraph, p::Float64, trials::Int)
    percolation_probabilities = zeros(Float64, trials)

    for trial in 1:trials
        # generate a percolated minimum spanning tree
        percolates = percolate_mst(graph, p)

        # record the percolation probability
        percolation_probabilities[trial] = percolates
    end

    return mean(percolation_probabilities)
end
