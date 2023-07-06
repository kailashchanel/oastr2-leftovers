#=
This file generates random points in the style of a GAMA csv file to be read by data_3d_distance.jl
=#


const NUM_GALAXIES::Int = 10000  # The number of galaxies to generate
const MAX_DISTANCE_MPC::Float64 = 1000  # The maximum distance galaxies can be, in megaparsecs

# Parameters for the G15 region
const RA_MAX::Float64 = 223.5
const RA_MIN::Float64 = 211.5

const DEC_MAX::Float64 = 3
const DEC_MIN::Float64 = -2


function generate_line()
    # The 0.0001 means that it generates a random number with 4 decimal places
    rand_ra = rand(RA_MIN:0.0001:RA_MAX)
    rand_dec = rand(DEC_MIN:0.0001:DEC_MAX)
    rand_dist = rand(0:0.0001:MAX_DISTANCE_MPC)

    return "$rand_ra,$rand_dec,$rand_dist\n"
end


function main()
    csv_file = "RA,DEC,dist\n"
    
    for _ in 1:NUM_GALAXIES
        csv_file *= generate_line()
    end

    file = open("rand_points.csv", "w")
    write(file, csv_file)
    close(file)
end


main()
