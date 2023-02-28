using Cosmology
using DataFrames
using CSV
using Plots

# Install packages
# using Pkg
# Pkg.add("Cosmology")
# Pkg.add("CSV")
# Pkg.add("Plots")
# Pkg.add("DataFrames")

Pkg.status()

#define cosmological model. For this example I will use the Planck 2015 
#cosmological parameters but this can be easily modified. 
c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

age(c, 0)
age(c, 1.42) #z = 1.42

data1 = CSV.read("path/to/data/file.csv", DataFrame)
data1

data1.age = age.(Ref(c), data1.Z_HELIO)

data1[!,"age"]

G02 = data1[((data1[!, "RA"] .< 38.8) .& (data1[!, "RA"].>30.2) .& (data1[!, "DEC"].<-3.72) .& (data1[!, "DEC"].>-10.25)),:]
G09 = data1[((data1[!, "RA"].<141.0) .& (data1[!, "RA"].>129.0) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
G12 = data1[((data1[!, "RA"].<186.0) .& (data1[!, "RA"].>174.0) .& (data1[!, "DEC"].<2.0) .& (data1[!, "DEC"].>-3.0)),:]
G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
G23 = data1[((data1[!, "RA"].<351.9) .& (data1[!, "RA"].>338.1) .& (data1[!, "DEC"].<-30.0) .& (data1[!, "DEC"].> -35.0)),:]

scatter(G02[!, "age"], G02[!, "RA"], G02[!, "DEC"])

# This saves the generated graph as a .png file in the current folder.
savefig("lookback-graph.png")