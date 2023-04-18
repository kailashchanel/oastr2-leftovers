#=
Generate an MST of the dataset. Note that this file generates a subset of the "actual" MST.
=#

import Pkg
Pkg.add(["Graphs", "SimpleGraphs", "GraphPlot", "CSV", "DataFrames", "Cosmology", "Plots", "SGtSNEpi", "GLMakie"])
using Graphs, SimpleGraphs, GraphPlot, CSV, DataFrames, Cosmology, Plots, GLMakie, SGtSNEpi

c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

data1 = CSV.read("GAMA_CZ5Unj.csv", DataFrame)
data1.age = age.(Ref(c), data1.Z_HELIO)

G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]

G15_matrix=Matrix(G15)

df = select(G15,([:RA,:DEC]))
    matrix_df=Matrix(df)
    t = matrix_df'

dfT = DataFrame(t, :auto)
matrix = Matrix(dfT)

@time g, dists = euclidean_graph(matrix[:,1:10000])

edgelist = collect(edges(g))

[[src(e), dst(e)] for e in edgelist]

@time k,dists = euclidean_graph(matrix[:,1:1000])

GLMakie.activate!()

y = sgtsnepi(k);
show_embedding(y;
  A = adjacency_matrix(g),        # show edges on embedding
  mrk_size = 1,                   # control node sizes
  lwd_in = 0.01, lwd_out = 0.001, # control edge widths
  edge_alpha = 0.4 )              # control edge transparency
