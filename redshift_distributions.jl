print("Install packages? (y/n) | ")
res = readline()

if cmp(res, "y") == 0
    using Pkg
    Pkg.add(["Distributions", "DataFrames", "Plots", "StatsPlots", "CSV"])
end

using DataFrames
using Distributions
using Plots
using StatsPlots
using CSV

data1 = CSV.read("/path/to/file/GAMA_CZ5Unj.csv", DataFrame)

G02 = data1[((data1[!, "RA"] .< 38.8) .& (data1[!, "RA"].>30.2) .& (data1[!, "DEC"].<-3.72) .& (data1[!, "DEC"].>-10.25)),:]
G09 = data1[((data1[!, "RA"].<141.0) .& (data1[!, "RA"].>129.0) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
G12 = data1[((data1[!, "RA"].<186.0) .& (data1[!, "RA"].>174.0) .& (data1[!, "DEC"].<2.0) .& (data1[!, "DEC"].>-3.0)),:]
G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
G23 = data1[((data1[!, "RA"].<351.9) .& (data1[!, "RA"].>338.1) .& (data1[!, "DEC"].<-30.0) .& (data1[!, "DEC"].> -35.0)),:]

data_arr = [G02, G09, G12, G15, G23]
name_arr = ["G02", "G09", "G12", "G15", "G23"]

for n in eachindex(data_arr)
    h = histogram(data_arr[n].Z_HELIO, linecolor = "blue", label="Data", color=:blue, normalize = true, title = string(name_arr[n], " Redshift Distribution"), xlabel = "Redshift", ylabel = "Frequency")

    dest = string("distributions/", name_arr[n],"-redshift-distribution.png")
    savefig(h, dest)

    println(string("\n", "Saved to '", dest, "'."))
end