using Clustering, Cosmology, DataFrames, CSV, Unitful, PlotlyJS

"Add dist values to the dataset"
function add_dist!(data)
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

    data.dist = comoving_radial_dist.(Ref(c), data.Z_HELIO)
end

function generate_graph(data, outputfilepath, name)
    raw_dist_arr = Vector{Float64}()
    data.raw_dist = Vector{Float64}(undef, length(data[!, "dist"]))
    
    for n in eachindex(data[!, "dist"])
        append!(raw_dist_arr, Unitful.ustrip(data[!, "dist"][n]))
    end
    
    data.raw_dist = raw_dist_arr

    data_matrix = Matrix{Float64}(undef, 3, 0)
    for n in eachindex(data[!, "RA"])
        A = [data[!, "RA"][n]; data[!, "DEC"][n]; data[!, "raw_dist"][n];]
        data_matrix = hcat(data_matrix, A)
    end

    println(size(data_matrix))

    R = kmeans(data_matrix, 10; maxiter=200, display=:iter)

    println("Centers of clusters:")
    println(R.centers)
    println("Number in each cluster:")
    println(counts(R))

    data.assignments = Vector{Float64}(undef, length(assignments(R)))
    data.assignments = assignments(R)

    p = plot(
        data,
        x=:raw_dist,
        y=:RA, 
        z=:DEC,
        color=:assignments,
        type="scatter3d", 
        mode="markers",
        Layout(
            title=string("RA, Dec, Radial Distance: ", name), 
            # xaxis=attr(title="Radial Distance"),
            # yaxis=attr(title="RA"),
            # zaxis=attr(title="Dec"),
            font=attr(
                family="Courier New",
                size=18,
            )
        )
    )
    
    open(string("./output/", outputfilepath), "w") do io
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

    # data1.age = age.(Ref(c), data1.Z_HELIO)

    println("Calculating radial distance...")
    add_dist!(data1)

    println("Calculating data groups...")
    G02 = data1[((data1[!, "RA"].<38.8) .& (data1[!, "RA"].>30.2) .& (data1[!, "DEC"].<-3.72) .& (data1[!, "DEC"].>-10.25)),:]
    G09 = data1[((data1[!, "RA"].<141.0) .& (data1[!, "RA"].>129.0) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
    G12 = data1[((data1[!, "RA"].<186.0) .& (data1[!, "RA"].>174.0) .& (data1[!, "DEC"].<2.0) .& (data1[!, "DEC"].>-3.0)),:]
    G15 = data1[((data1[!, "RA"].<223.5) .& (data1[!, "RA"].>211.5) .& (data1[!, "DEC"].<3.0) .& (data1[!, "DEC"].>-2.0)),:]
    G23 = data1[((data1[!, "RA"].<351.9) .& (data1[!, "RA"].>338.1) .& (data1[!, "DEC"].<-30.0) .& (data1[!, "DEC"].> -35.0)),:]

    println("Generating graphs...")
    generate_graph(G02, "g02-clustered-10.html", "G02")
    # generate_graph(G09, "g09-radial.html", "G09")
    # generate_graph(G12, "g12-radial.html", "G12")
    # generate_graph(G15, "g15-radial.html", "G15")
    # generate_graph(G23, "g23-radial.html", "G23")

    println("Program complete.")
end


main()