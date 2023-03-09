import Pkg
Pkg.add(["Graphs", "SimpleWeightedGraphs", "GraphPlot"])

using Graphs, SimpleWeightedGraphs, GraphPlot


graph = SimpleWeightedGraph(1)


function add_edges(filepath)
    lines = countlines(filepath)
    i = 0
    for row::String in eachline(filepath)
        i += 1
        if row[1] == 'v'  # ignore csv header
            continue
        end

        splitted = split(row, ",")
        
        v1::Int = parse(Int, splitted[1])
        v2::Int = parse(Int, splitted[2])
        dist::Float32 = parse(Float32, splitted[3])
        
        if v2 > nv(graph) 
            add_vertex!(graph)
        end

        add_edge!(graph, v1, v2, dist)
        if i % 10000 == 0
            println(string(i / lines * 100, "% done"))
        end
    end
end


function main()
    for file in readdir("G02")
        println("Calculating file $file")
        add_edges("G02/$file")
        println("Added eges for file $file")
    end

    savegraph("G02.lgz", graph)
end


@time main()
