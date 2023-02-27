using Pkg

Pkg.status()

using Random

X = rand(100)
Y = rand(100)
Z = rand(100)

scatter3D(X,Y,Z)

using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")

using CSV
using DataFrames


data1 = CSV.read("/Users/sophiestrey/Desktop/data/ASTROPHYSICS/leftovers/GAMA_CZ5Unj.csv", DataFrame)

x = data1[!, "RA"]
y = data1[!, "DEC"]
scatter(x,y)

G02 = data1[((data1[!, "RA"] .< 38.8) .& (data1[!, "RA"].>30.2) .& (data1[!, "DEC"].<-3.72) .& (data1[!, "DEC"].>-10.25)),:]


typeof(data1)

Pkg.add("Plots")
using Plots

scatter(G02[!,"RA"], G02[!, "DEC"])

scatter3D(G02[!, "Z_HELIO"], G02[!, "RA"], G02[!, "DEC"])

G09 = data1[((data1[!, "RA"].<141.0) .& (data1[!, "RA"].>129.0) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
G12 = data1[((data1[!, "RA"].<186.0) .& (data1[!, "RA"].>174.0) .& (data1[!, "DEC"].<2.0) .& (data1[!, "DEC"].>-3.0)),:]
G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
G23 = data1[((data1[!, "RA"].<351.9) .& (data1[!, "RA"].>338.1) .& (data1[!, "DEC"].<-30.0) .& (data1[!, "DEC"].> -35.0)),:]

scatter3D(G02[!, "Z_HELIO"], G02[!, "RA"], G02[!, "DEC"])
scatter3D(G09[!, "Z_HELIO"], G09[!, "RA"], G09[!, "DEC"])
scatter3D(G12[!, "Z_HELIO"], G12[!, "RA"], G12[!, "DEC"])
scatter3D(G15[!, "Z_HELIO"], G15[!, "RA"], G15[!, "DEC"])
scatter3D(G23[!, "Z_HELIO"], G23[!, "RA"], G23[!, "DEC"])

scatter3D(G02[!, "Z_HELIO"], G02[!, "RA"], G02[!, "DEC"])

scatter3D(G09[!, "Z_HELIO"], G09[!, "RA"], G09[!, "DEC"])

scatter3D(G12[!, "Z_HELIO"], G12[!, "RA"], G12[!, "DEC"])

scatter3D(G15[!, "Z_HELIO"], G15[!, "RA"], G15[!, "DEC"])

scatter3D(G23[!, "Z_HELIO"], G23[!, "RA"], G23[!, "DEC"])


