using Pkg
Pkg.add(["Cosmology", "DataFrames", "CSV", "Unitful", "PlotlyJS"])

using Cosmology, DataFrames, CSV, Unitful, PlotlyJS


function generate_graph(data, outputfilepath)
    raw_age_arr = Vector{Float64}()
    data.raw_age = Vector{Float64}(undef, length(data[!, "age"]))
    
    for n in eachindex(data[!, "age"])
        append!(raw_age_arr, Unitful.ustrip(data[!, "age"][n]))
    end
    
    data.raw_age = raw_age_arr
    
    p = plot(
        data,
        x=:raw_age,
        y=:RA, 
        z=:DEC,
        type="scatter3d", 
        mode="markers",
        Layout(
            title="RA & Dec Through Time", 
            xaxis_title="Lookback Time (Gyr)",
            yaxis_title="RA", 
            zaxis_title="Dec", 
            font=attr(
                family="Courier New",
                size=18,
            )
        )
    )
    
    open(outputfilepath, "w") do io
        PlotlyBase.to_html(io, p.plot)
    end
end


function main()
    #define cosmological model. For this example I will use the Planck 2015 
    #cosmological parameters but this can be easily modified. 

    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

    age(c, 0)
    age(c, 1.42) #z = 1.42

    data1 = CSV.read("GAMA_CZ5Unj.csv", DataFrame)
    data1

    data1.age = age.(Ref(c), data1.Z_HELIO)

    data1[!,"age"]

    G02 = data1[((data1[!, "RA"].<38.8) .& (data1[!, "RA"].>30.2) .& (data1[!, "DEC"].<-3.72) .& (data1[!, "DEC"].>-10.25)),:]
    G09 = data1[((data1[!, "RA"].<141.0) .& (data1[!, "RA"].>129.0) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
    G12 = data1[((data1[!, "RA"].<186.0) .& (data1[!, "RA"].>174.0) .& (data1[!, "DEC"].<2.0) .& (data1[!, "DEC"].>-3.0)),:]
    G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
    G23 = data1[((data1[!, "RA"].<351.9) .& (data1[!, "RA"].>338.1) .& (data1[!, "DEC"].<-30.0) .& (data1[!, "DEC"].> -35.0)),:]

    generate_graph(G02, "./g02.html")
    generate_graph(G09, "./g09.html")
    generate_graph(G12, "./g12.html")
    generate_graph(G15, "./g15.html")
    generate_graph(G23, "./g23.html")
end


main()
