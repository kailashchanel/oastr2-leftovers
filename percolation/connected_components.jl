using LightGraphs, GraphPlot, Random, Plots

#sample distance matrix
dist_matrix = [0 1 2 3 4; 
                1 0 1 2 3; 
                2 1 0 1 2; 
                3 2 1 0 1; 
                4 3 2 1 0]


g = Graph(convert(Matrix{Bool}, dist_matrix .<= 1.5))

#testing functions 
deg_cent = degree_centrality(g)
println("Degree centrality: ", deg_cent)

bet_cent = betweenness_centrality(g)
println("Betweenness centrality: ", bet_cent)

comms = connected_components(g)
println("Communities: ", comms)

GraphPlot.gplot(g)
