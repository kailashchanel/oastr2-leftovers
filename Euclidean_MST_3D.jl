import Pkg
Pkg.add("Graphs")
Pkg.add("SimpleGraphs")
Pkg.add("GraphPlot")

using Graphs
using SimpleGraphs
using GraphPlot

points = rand(3,10) # n=10 vertices in R^3 (3 dimesions)

euclidean_graph

pts = rand(3,10);
g, dists = euclidean_graph(pts, p=1, bc=:periodic);
g

GraphPlot.gplot(g)

dists

mst = Graphs.boruvka_mst(g)   

GraphPlot.gplot(mst)

using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")

using CSV
using DataFrames

Pkg.add("Cosmology")

using Cosmology
c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

data1 = CSV.read("/Users/sophiestrey/Desktop/data/ASTROPHYSICS/leftovers/GAMA_CZ5Unj.csv", DataFrame)
data1.age = age.(Ref(c), data1.Z_HELIO)
data1

G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]

pts = G15[G15[!, "RA"]], G15[G15[!, "DEC"]]

g, dists = euclidean_graph(pts, p=2, bc=:periodic);
g

G15_matrix=Matrix(G15)

df = select(G15,([:RA,:DEC]))
matrix_df=Matrix(df)

transpose(matrix_df)
t = matrix_df'

dfT = DataFrame(t, :auto)
matrix = Matrix(dfT)

g, dists = euclidean_graph(points);
g

@time g,dists = euclidean_graph(matrix)


