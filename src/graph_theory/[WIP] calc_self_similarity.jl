using FractalDimensions
using DynamicalSystemsBase
using DelimitedFiles

using Graphs, SimpleGraphs, GraphPlot, CSV, DataFrames, Cosmology, Plots, GLMakie, SGtSNEpi

c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

data1 = CSV.read("/Users/sophiestrey/Desktop/GAMA_CZ5Unj.csv", DataFrame)
data1.age = age.(Ref(c), data1.Z_HELIO)

G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]


G15_matrix=Matrix(G15)

df = select(G15,([:RA,:DEC,:Z_HELIO]))
    matrix_df=Matrix(df)
    t = matrix_df'

dfT = DataFrame(t, :auto)
matrix = Matrix(dfT)

#X = StateSpaceSet(matrix)
e = 2 .^ (-15:0.5:5) #semi-random guess
Cs = correlationsum(matrix, e; show_progress = false)

