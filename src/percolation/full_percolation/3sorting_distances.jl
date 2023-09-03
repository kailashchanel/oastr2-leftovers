#in use for analysis on complete graph. Do NOT use for MST. 
#creates struct 
using JLD2, Plots

struct GalaxyDistance
    idx1::Int32
    idx2::Int32
    distance::Float32
end

distance_list = load_object("distance_list_old.jld2")

Base.isless(a::GalaxyDistance, b::GalaxyDistance) = a.distance < b.distance
Base.isequal(a::GalaxyDistance, b::GalaxyDistance) = a.distance == b.distance

function sort_distances!(distances)
    sort!(distances, by = x -> x.distance)
end

sort_distances!(distance_list_young)

save_object("sorted_distance_list_young.jld2", distance_list_young)

histogram([x.distance for x in distance_list],bins=1000,legend=false,xlabel="Distance (Mpc)",ylabel="Number of Galaxies",title="Distance Histogram")
savefig("distance_histogram.png")
