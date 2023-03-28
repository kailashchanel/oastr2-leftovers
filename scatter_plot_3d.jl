print("Install packages? (y/n) | ")
res = readline()

if cmp(res, "y") == 0
    using Pkg
    Pkg.add(["Cosmology", "DataFrames", "CSV", "Unitful", "PlotlyJS"])
end

using Cosmology
using DataFrames
using CSV
using Unitful
using PlotlyJS

#define cosmological model. For this example I will use the Planck 2015 
#cosmological parameters but this can be easily modified. 

c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

age(c, 0)
age(c, 1.42) #z = 1.42

data1 = CSV.read("C:\\Users\\kaila\\Downloads\\GAMA_CZ5Unj.csv", DataFrame)
data1

data1.age = age.(Ref(c), data1.Z_HELIO)

data1[!,"age"]

G02 = data1[((data1[!, "RA"] .< 38.8) .& (data1[!, "RA"].>30.2) .& (data1[!, "DEC"].<-3.72) .& (data1[!, "DEC"].>-10.25)),:]
G09 = data1[((data1[!, "RA"].<141.0) .& (data1[!, "RA"].>129.0) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
G12 = data1[((data1[!, "RA"].<186.0) .& (data1[!, "RA"].>174.0) .& (data1[!, "DEC"].<2.0) .& (data1[!, "DEC"].>-3.0)),:]
G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
G23 = data1[((data1[!, "RA"].<351.9) .& (data1[!, "RA"].>338.1) .& (data1[!, "DEC"].<-30.0) .& (data1[!, "DEC"].> -35.0)),:]

raw_age_arr = Vector{Float64}()
G02.raw_age = Vector{Float64}(undef, length(G02[!, "age"]))

for n in eachindex(G02[!, "age"])
    append!(raw_age_arr, Unitful.ustrip(G02[!, "age"][n]))
end

G02.raw_age = raw_age_arr

p = plot(
    G02,
    x=:raw_age,
    y=:RA, 
    z=:DEC,
    type="scatter3d", 
    mode="markers",
    Layout(
        title="RA & Dec Through Time", 
        xaxis_title="RA",
        yaxis_title="Lookback Time (Gyr)", 
        zaxis_title="Declination", 
        font=attr(
            family="Courier New",
            size=18,
        )
    )
)

open("./output/scatter3d.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end