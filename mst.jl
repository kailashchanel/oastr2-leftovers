import Pkg
Pkg.add("Graphs")
Pkg.add("GraphPlot")

import Graphs
import GraphPlot

g = Graphs.ladder_graph(6)

mst = Graphs.boruvka_mst(g)

new_graph = Graphs.cycle_graph(0)
Graphs.add_vertices!(new_graph, Graphs.nv(g))

for edge in mst.mst
    Graphs.add_edge!(new_graph, edge)
end

GraphPlot.gplot(new_graph)
