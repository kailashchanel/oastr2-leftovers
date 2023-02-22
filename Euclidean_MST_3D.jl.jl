import Pkg
Pkg.add("Graphs")
Pkg.add("SimpleGraphs")
Pkg.add("GraphPlot")
Pkg.add(Pkg.PackageSpec(name="Graphs",rev="master"))

import Graphs
import SimpleGraphs
import GraphPlot

points = rand(3,10) # n=10 vertices in R^3 (3 dimesions)

g, dists = euclidean_graph(points, p=2, bc=:periodic) 
g

min = boruvka_mst(g)

for edge in min.mst
    Graphs.add_edge!(new_graph, edge)
end

min = boruvka_mst(g)
GraphPlot.gplot(new_graph)


