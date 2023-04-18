using CSV, DataFrames, Cosmology, JLD2

function spherical_to_cart(ra, dec, dist)
    # Convert spherical coordinates to cartesian coordinates
    # RA and Dec are in degrees
    # dist is in Mpc
    # Returns x, y, z in Mpc
    x = dist * cos(deg2rad(dec)) * cos(deg2rad(ra))
    y = dist * cos(deg2rad(dec)) * sin(deg2rad(ra))
    z = dist * sin(deg2rad(dec))
    return x, y, z
end

data = CSV.read("GAMA_CZ5Unj.csv", DataFrame)

println("Calculating radial distance")
function add_dist!(data)
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

    data.dist = comoving_radial_dist.(Ref(c), data.Z_HELIO)
end
add_dist!(data)

println("Calculating cosmological age")
function add_age!(data)
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

    data.age = age.(Ref(c), data.Z_HELIO)
end
add_age!(data)

println("Calculating data groups")
# G02 = data[((data[!, "RA"] .< 38.8) .& (data[!, "RA"].>30.2) .& (data[!, "DEC"].<-3.72) .& (data[!, "DEC"].>-10.25)),:]
# G09 = data[((data[!, "RA"].<141.0) .& (data[!, "RA"].>129.0) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
# G12 = data[((data[!, "RA"].<186.0) .& (data[!, "RA"].>174.0) .& (data[!, "DEC"].<2.0) .& (data[!, "DEC"].>-3.0)),:]
G15 = data[((data[!, "RA"].<223.5) .& (data[!, "RA"].>211.5) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
# G23 = data[((data[!, "RA"].<351.9) .& (data[!, "RA"].>338.1) .& (data[!, "DEC"].<-30.0) .& (data[!, "DEC"].> -35.0)),:]

ra = G15[1:4:end, "RA"]
dec = G15[1:4:end, "DEC"]
distunit = G15[1:4:end, "dist"]
ageunit = G15[1:4:end, "age"]

#remove units
dist = Real.(distunit ./ oneunit.(distunit))
cage = Real.(ageunit ./ oneunit.(ageunit))
datalen = length(ra)

# Create a struct to store the distance data
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
        # RA and Dec are in degrees, dist is in Mpc
        x1, y1, z1 = spherical_to_cart(ra[i], dec[i], dist[i])
        x2, y2, z2 = spherical_to_cart(ra[j], dec[j], dist[j])
        distance = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
        
        push!(distance_list,GalaxyDistance(i,j,distance))
    end
end

save_object("distance_list025.jld2", distance_list)

histogram([d.distance for d in distance_list], bins=100, xlabel="Distance (Mpc)", ylabel="Count")
