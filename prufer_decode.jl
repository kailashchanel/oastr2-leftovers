# !! WARNING: THIS CODE DOES NOT WORK !!
# I tried making it work and failed

import Pkg
Pkg.add("Graphs")
Pkg.add("GraphPlot")

import Graphs
import GraphPlot


struct EnhancedNumber
    value::Int64
    available::Bool
end


# I used this video as a resource for this algorithm:
# https://www.youtube.com/watch?v=p-g9ZV3xuGc&ab_channel=ShannonDietz
function decode(l, s)
    sort!(l)
    sort!(s)

    graph = Graphs.path_graph(0)  # not sure of any other way to init an empty graph
    
    node_dict = Dict()

    for i in eachindex(s)
        if !haskey(node_dict, s[i].value)
            Graphs.add_vertex!(graph)
            node_dict[Graphs.nv(graph)] = i
        end
        
        for j in eachindex(l)
            if l[j].value < s[i].value || !l[j].available
                continue
            end
            
            if !haskey(node_dict, l[j])
                Graphs.add_vertex!(graph)
                node_dict[Graphs.nv(graph)] = i
            end

            Graphs.add_edge!(graph, s[i].value, l[j].value)
            s[i].available = false
            l[j].available = false
            break
        end
    end

    return graph
end

graph = decode(
    [
        EnhancedNumber(1, true),
        EnhancedNumber(2, true),
        EnhancedNumber(3, true),
        EnhancedNumber(4, true),
        EnhancedNumber(5, true),
        EnhancedNumber(6, true),
        EnhancedNumber(7, true)
    ],
    [
        EnhancedNumber(1, true),
        EnhancedNumber(1, true),
        EnhancedNumber(3, true),
        EnhancedNumber(5, true),
        EnhancedNumber(5, true)
    ]
)
GraphPlot.gplot(graph)
