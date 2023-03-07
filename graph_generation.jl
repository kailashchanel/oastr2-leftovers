import Pkg
Pkg.add("Graphs")
Pkg.add("SimpleGraphs")
Pkg.add("GraphPlot")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Cosmology")

using Graphs
using SimpleGraphs
using GraphPlot
using CSV
using DataFrames
using Cosmology

c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

data1 = CSV.read("GAMA_CZ5Unj.csv", DataFrame)
data1.age = age.(Ref(c), data1.Z_HELIO)
data1

G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]

G15_matrix=Matrix(G15)

df = select(G15,([:RA,:DEC]))
    matrix_df=Matrix(df)
    t = matrix_df'

dfT = DataFrame(t, :auto)
matrix = Matrix(dfT)

@time g,dists = euclidean_graph(matrix[:,1:10000])

g

Pkg.add("Plots")
using Plots

edgelist = collect(edges(g))

[[src(e), dst(e)] for e in edgelist]

@time k,dists = euclidean_graph(matrix[:,1:1000])

Pkg.add("SGtSNEpi")
Pkg.add("GLMakie")

using GLMakie, SGtSNEpi

GLMakie.activate!()

y = sgtsnepi(k);
show_embedding(y;
  A = adjacency_matrix(g),        # show edges on embedding
  mrk_size = 1,                   # control node sizes
  lwd_in = 0.01, lwd_out = 0.001, # control edge widths
  edge_alpha = 0.4 )             # control edge transparency


