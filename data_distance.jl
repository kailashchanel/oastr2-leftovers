#=
A multithreaded, brute-force program to calculate the distance matrix.
Calculating distance matrices will take less than 1 hour with 20 threads at
4.8 GHz.
=#

using Pkg
Pkg.add(["CSV", "DataFrames", "Cosmology", "Unitful"])
using CSV, DataFrames, Cosmology, Unitful


function write_distances(distance_data, filestream)
    for row in distance_data
        write(filestream, string(Int(row[1]), ",", Int(row[2]), ",", row[3], "\n"))
    end
end


"Calculate the distance from one point and every other point. Will write the output to filename"
function calculate_distances(data, filename, startindex, endindex)
    println("Start function with filename $filename startindex $startindex endindex $endindex")
    distances = []

    filestream = open(filename, "a")
    write(filestream, "v1,v2,dist\n")  # set up csv header

    for i in range(startindex, endindex)
        for j in range(i + 1, size(data)[1] - 1)
            distance = sqrt(
                (data[i, "RA"] - data[j, "RA"])^2 + 
                (data[i, "DEC"] - data[j, "DEC"])^2 +
                (ustrip(data[i, "age"]) - ustrip(data[j, "age"]))^2
            )
            
            push!(distances, [i, j, distance])            
        end

        write_distances(distances, filestream)
        distances = []
    end

    write_distances(distances, filestream)
    println(string("Done with range ", startindex, " to ", endindex))
    close(filestream)
end


"Add age values to the dataset"
function add_age!(data)
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

    data.age = age.(Ref(c), data.Z_HELIO)
end


"Does the same thing as the function calculate_distnaces however it is multithreaded and writes to multiple different output files"
function calculate_distances_multithreaded(data)
    startindex = 1
    jump = size(data)[1] / (Threads.nthreads() + 1)
    endindex = jump

    for i in 1:Threads.nthreads()
        println("Spawning thread with startindex $startindex and endindex $endindex, i = $i")
        @Threads.spawn calculate_distances(data, string("output", i, ".csv"), Int(ceil(startindex)), Int(ceil(endindex)))

        sleep(0.1)

        startindex += jump
        endindex += jump
    end
end


function main()
    println("Loading data")
    data = CSV.read("GAMA_CZ5Unj.csv", DataFrame)

    println("Calculating age")
    add_age!(data)

    println("Calculating data groups")
    G02 = data[((data[!, "RA"] .< 38.8) .& (data[!, "RA"].>30.2) .& (data[!, "DEC"].<-3.72) .& (data[!, "DEC"].>-10.25)),:]
    G09 = data[((data[!, "RA"].<141.0) .& (data[!, "RA"].>129.0) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
    G12 = data[((data[!, "RA"].<186.0) .& (data[!, "RA"].>174.0) .& (data[!, "DEC"].<2.0) .& (data[!, "DEC"].>-3.0)),:]
    G15 = data[((data[!, "RA"].<223.5) .& (data[!, "RA"].>211.5) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
    G23 = data[((data[!, "RA"].<351.9) .& (data[!, "RA"].>338.1) .& (data[!, "DEC"].<-30.0) .& (data[!, "DEC"].> -35.0)),:]

    println("Calculating distance matrix")

    calculate_distances_multithreaded(G02)
    
    println("Done")
end


main()
