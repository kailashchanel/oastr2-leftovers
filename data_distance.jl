# watched this video on dataframes
# https://www.youtube.com/watch?v=Pt8Iz4Udg2E&ab_channel=TheJuliaProgrammingLanguage

using Pkg
Pkg.add(["CSV", "DataFrames", "Cosmology", "Unitful"])
using CSV, DataFrames, Cosmology, Unitful


function make_distance_df(data)
    distance_df = DataFrame(index1=[], index2=[], distance=[])

    for i in range(1, size(data)[1])
        for j in range(i + 1, size(data)[1] - 1)
            distance = sqrt(
                (data[i, "RA"] - data[j, "RA"])^2 + 
                (data[i, "DEC"] - data[j, "DEC"])^2 +
                ustrip((data[i, "age"] - (data[j, "age"])^2))
            )

            push!(distance_df, [i, j, distance])
        end

        println(string(i / size(data)[1], "% complete"))
    end

    return distance_df
end


function add_age!(data)
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

    data.age = age.(Ref(c), data.Z_HELIO)
end


function main()
    println("Loading data")
    data = CSV.read("GAMA_CZ5Unj.csv", DataFrame)

    println("Calculating age")
    add_age!(data)

    println("Calculating data groups")
    G09 = data[((data[!, "RA"].<141.0) .& (data[!, "RA"].>129.0) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
    G12 = data[((data[!, "RA"].<186.0) .& (data[!, "RA"].>174.0) .& (data[!, "DEC"].<2.0) .& (data[!, "DEC"].>-3.0)),:]
    G15 = data[((data[!, "RA"].<223.5) .& (data[!, "RA"].>211.5) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
    G23 = data[((data[!, "RA"].<351.9) .& (data[!, "RA"].>338.1) .& (data[!, "DEC"].<-30.0) .& (data[!, "DEC"].> -35.0)),:]

    println("Calculating distance matrix")
    distance_df = make_distance_df(G09)
    println("Done")

    println(distance_df)
end

main()
