using CSV, DataFrames, Cosmology, JLD2, Plots

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
    c_age = age.(Ref(c), data.Z_HELIO)
    # add age without units
    data.age = c_age ./ oneunit.(c_age)
end
add_age!(data)

println("Calculating data groups")
# G02 = data[((data[!, "RA"] .< 38.8) .& (data[!, "RA"].>30.2) .& (data[!, "DEC"].<-3.72) .& (data[!, "DEC"].>-10.25)),:]
# G09 = data[((data[!, "RA"].<141.0) .& (data[!, "RA"].>129.0) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
# G12 = data[((data[!, "RA"].<186.0) .& (data[!, "RA"].>174.0) .& (data[!, "DEC"].<2.0) .& (data[!, "DEC"].>-3.0)),:]
G15_unfiltered = data[((data[!, "RA"].<223.5) .& (data[!, "RA"].>211.5) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0) .& (data[!,"age"].>7.5)),:]
# G23 = data[((data[!, "RA"].<351.9) .& (data[!, "RA"].>338.1) .& (data[!, "DEC"].<-30.0) .& (data[!, "DEC"].> -35.0)),:]

G15 = filter(row -> row.age != 0, G15_unfiltered)

G15_old = G15[(G15[!,"age"].>11,:) & G15[!, "age"].<13]
G15_young = G15[(G15[!,"age"].<11,:) & G15[!, "age"].<13]
G15_middle = G15[(G15[!, "age"])] #fill in ranges

ra_old = G15_old[:, "RA"]
dec_old = G15_old[:, "DEC"]
distunit_old = G15_old[:, "dist"]
c_age_old = G15_old[:, "age"]

ra_young = G15_young[:, "RA"]
dec_young = G15_young[:, "DEC"]
distunit_young = G15_young[:, "dist"]
c_age_young = G15_young[:, "age"]

ra_middle = G15_middle[:, "RA"]
dec_middle = G15_middle[:, "DEC"]
distunit_middle = G15_middle[:, "dist"]
c_age_middle = G15_middle[:, "age"]

#remove units
dist_old = Real.(distunit_old ./ oneunit.(distunit_old))
dist_young = Real.(distunit_young ./ oneunit.(distunit_young))
dist_middle = Real.(distunit_middle ./ oneunit.(distunit_middle))
datalen_old = length(ra_old)
datalen_young = length(ra_young)
datalen_middle = length(ra_middle)

# Plot the age histogram
histogram(c_age_young, bins=500, xlabel="Cosmological age (Gyr)", ylabel="Count",legend = false)
savefig("young_age_hist.png")
histogram(c_age_middle, bins=500, xlabel="Cosmological age (Gyr)", ylabel="Count",legend = false)
savefig("middle_age_hist.png")
histogram(c_age_old, bins=500, xlabel="Cosmological age (Gyr)", ylabel="Count",legend = false)
savefig("old_age_hist.png")

# Create a struct to store the distance data
struct GalaxyDistance
    idx1::Int32
    idx2::Int32
    distance::Float32
end

distance_list_old = GalaxyDistance[]

for i in 1:datalen_old-1
    println("i = $i")
    for j in i+1:datalen_old
        # Find the separation in the sky
        # RA and Dec are in degrees, dist is in Mpc
        x1, y1, z1 = spherical_to_cart(ra_old[i], dec_old[i], dist_old[i])
        x2, y2, z2 = spherical_to_cart(ra_old[j], dec_old[j], dist_old[j])
        distance = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)

        push!(distance_list_old,GalaxyDistance(i,j,distance))
    end
end

distance_list_young = GalaxyDistance[]

for i in 1:datalen_young-1
    println("i = $i")
    for j in i+1:datalen_young
        # Find the separation in the sky
        # RA and Dec are in degrees, dist is in Mpc
        x1, y1, z1 = spherical_to_cart(ra_young[i], dec_young[i], dist_young[i])
        x2, y2, z2 = spherical_to_cart(ra_young[j], dec_young[j], dist_young[j])
        distance = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)

        push!(distance_list_young,GalaxyDistance(i,j,distance))
    end
end

distance_list_middle = GalaxyDistance[]

for i in 1:datalen_middle-1
    println("i = $i")
    for j in i+1:datalen_middle
        # Find the separation in the sky
        # RA and Dec are in degrees, dist is in Mpc
        x1, y1, z1 = spherical_to_cart(ra_middle[i], dec_middle[i], dist_middle[i])
        x2, y2, z2 = spherical_to_cart(ra_middle[j], dec_middle[j], dist_middle[j])
        distance = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)

        push!(distance_list_middle, GalaxyDistance(i,j,distance))
    end
end

save_object("distance_list_old.jld2", distance_list_old)
save_object("distance_list_young.jld2", distance_list_young)
save_object("distance_list_middle.jld2", distance_list_middle)

histogram([d.distance for d in distance_list_old], bins=100, xlabel="Distance old galaxies (Mpc)", ylabel="Count",legend = false)
savefig("distance_list_old.png")
histogram([d.distance for d in distance_list_young], bins=100, xlabel="Distance young galaxies (Mpc)", ylabel="Count",legend = false)
savefig("distance_list_young.png")
histogram([d.distance for d in distance_list_middle], bins=100, xlabel="Distance middle galaxies (Mpc)", ylabel="Count",legend = false)
savefig("distance_list_middle.png")