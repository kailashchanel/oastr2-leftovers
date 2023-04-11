using JLD2

struct GalaxyDistance
    idx1::Int32
    idx2::Int32
    distance::Float32
end

distances = load_object("distance_list.jld2")
