using CSV, DataFrames, Cosmology, JLD2

data = CSV.read("GAMA_CZ5Unj.csv", DataFrame)

println("Calculating radial distance")
function add_dist!(data)
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

    data.dist = comoving_radial_dist.(Ref(c), data.Z_HELIO)
end
add_dist!(data)

println("Calculating data groups")
# G02 = data[((data[!, "RA"] .< 38.8) .& (data[!, "RA"].>30.2) .& (data[!, "DEC"].<-3.72) .& (data[!, "DEC"].>-10.25)),:]
# G09 = data[((data[!, "RA"].<141.0) .& (data[!, "RA"].>129.0) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
# G12 = data[((data[!, "RA"].<186.0) .& (data[!, "RA"].>174.0) .& (data[!, "DEC"].<2.0) .& (data[!, "DEC"].>-3.0)),:]
G15 = data[((data[!, "RA"].<223.5) .& (data[!, "RA"].>211.5) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
# G23 = data[((data[!, "RA"].<351.9) .& (data[!, "RA"].>338.1) .& (data[!, "DEC"].<-30.0) .& (data[!, "DEC"].> -35.0)),:]

ra = G15[:, "RA"]
dec = G15[:, "DEC"]
distunit = G15[:, "dist"]
dist = Real.(distunit ./ oneunit.(distunit))
datalen = length(ra)

struct GalaxyDistance
    idx1::Int32
    idx2::Int32
    distance::Float32
end

distance_list = GalaxyDistance[]

for i in 1:datalen-1
    println("i = $i")
    for j in i+1:datalen
        # Find the separation in the sky
        # RA and Dec are in degrees
        deltaRA = abs(ra[i] - ra[j])
        deltaDec = abs(dec[i] - dec[j])
        skySep = sqrt(deltaRA^2 + deltaDec^2)

        # Find transverse separation. The two distances are averaged because the formula requires one distance.
        # The two distances averaged is a compromise.
        sep_transverse = 2 * sin(deg2rad(skySep / 2)) * ((dist[i] + dist[j]) / 2)
        
        distance = sqrt(sep_transverse^2 + (abs(dist[i] - abs(dist[j]))^2))

        push!(distance_list,GalaxyDistance(i,j,distance_matrix[i,j]))
    end
end

save_object("distance_list.jld2", distance_list)

