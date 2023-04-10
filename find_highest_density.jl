import Pkg
Pkg.add(["CSV", "Cosmology", "DataFrames", "Unitful"])
using CSV, Cosmology, DataFrames, Unitful


struct area3
    x_min
    x_max
    y_min
    y_max
    z_min
    z_max
end


function add_xyz(data)
    xyz::DataFrame = DataFrame(x=[], y=[], z=[])

    for datapoint in eachrow(data)
        ra_hyp = cos(deg2rad(datapoint.DEC)) * datapoint.dist
        x = cos(deg2rad(datapoint.RA)) * ra_hyp
        y = sin(deg2rad(datapoint.DEC))
        z = sin(deg2rad(datapoint.RA)) * datapoint.dist

        push!(xyz, [x, y, z])
    end

    return hcat(data, xyz)
end


"Add dist values to the dataset"
function add_dist!(data)
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

    data.dist = comoving_radial_dist.(Ref(c), data.Z_HELIO)
end

function add_axes!(data)
    # Find the separation in the sky
    # RA and Dec are in degrees
    deltaRA = abs(data[i, "RA"] - data[j, "RA"])
    deltaDec = abs(data[i, "DEC"] - data[j, "DEC"])
    skySep = sqrt(deltaRA^2 + deltaDec^2)

    # Find transverse separation. The two distances are averaged because the formula requires one distance.
    # The two distances averaged is a compromise.
    sep_transverse = 2 * sin(deg2rad(skySep / 2)) * ((data[i, "dist"] + data[j, "dist"]) / 2)
    
    dist = sqrt(sep_transverse^2 + (abs(data[i, "dist"] - abs(data[j, "dist"]))^2))

    push!(distances, [i, j, dist])
end


"Returns units of MPc^3/galaxies"
function calculate_density(data::DataFrame, x_min, x_max, y_min, y_max, z_min, z_max)
    area = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)

    galaxiesWithinArea::Int = 0

    x_min = ustrip(x_min)
    x_max = ustrip(x_max)
    y_min = ustrip(y_min)
    y_max = ustrip(y_max)
    z_min = ustrip(z_min)
    z_max = ustrip(z_max)

    for row in eachrow(data)
        if ustrip(row.x) >= x_min && ustrip(row.x) <= x_max &&
            ustrip(row.y) >= y_min && ustrip(row.y) <= y_max &&
            ustrip(row.z) >= z_min && ustrip(row.z) <= z_max
            
            galaxiesWithinArea += 1
        end
    end

    return area / galaxiesWithinArea
end


"finalArea parameter is the final area in megaparsecs.
Having this value to be lower is more accurate but more computationally expensive. An infinite loop occurs when
finalArea is 0 or less."
function findHighestDensityArea(data, finalArea)
    x_min = findmin(data.x)[1]
    x_max = findmax(data.x)[1]
    y_min = findmin(data.y)[1]
    y_max = findmax(data.y)[1]
    z_min = findmin(data.z)[1]
    z_max = findmax(data.z)[1]

    # From here on this function is spaghetti code but I'm not sure how to clean it up
    x_lower = x_min
    x_upper = x_max

    while (ustrip(x_upper - x_lower) > finalArea)
        middle = (x_upper + x_lower) / 2

        left = calculate_density(
            data, x_lower, middle, y_min, y_max, z_min, z_max
        )

        right = calculate_density(
            data, middle, x_upper, y_min, y_max, z_min, z_max
        )

        if left > right
            x_upper = middle
        else
            x_lower = middle
        end
    end

    y_lower = y_min
    y_upper = y_max

    while (ustrip(y_upper - y_lower) > finalArea)
        middle = (y_upper + y_lower) / 2

        left = calculate_density(
            data, y_lower, middle, y_min, y_max, z_min, z_max
        )

        right = calculate_density(
            data, middle, y_upper, y_min, y_max, z_min, z_max
        )

        if left > right
            y_upper = middle
        else
            y_lower = middle
        end
    end

    z_lower = z_min
    z_upper = z_max

    while (ustrip(z_upper - z_lower) > finalArea)
        middle = (z_upper + z_lower) / 2

        left = calculate_density(
            data, z_lower, middle, z_min, z_max, z_min, z_max
        )

        right = calculate_density(
            data, middle, z_upper, z_min, z_max, z_min, z_max
        )

        if left > right
            z_upper = middle
        else
            z_lower = middle
        end
    end

    return area3(x_lower, x_upper, y_lower, y_upper, z_lower, z_upper)
end


function main()
    println("Loading data")
    data = CSV.read("GAMA_CZ5Unj.csv", DataFrame)

    println("Calculating radial distance")
    add_dist!(data)
    data = add_xyz(data)


    println("Calculating data groups")
    G02 = data[((data[!, "RA"] .< 38.8) .& (data[!, "RA"].>30.2) .& (data[!, "DEC"].<-3.72) .& (data[!, "DEC"].>-10.25)),:]
    G09 = data[((data[!, "RA"].<141.0) .& (data[!, "RA"].>129.0) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
    G12 = data[((data[!, "RA"].<186.0) .& (data[!, "RA"].>174.0) .& (data[!, "DEC"].<2.0) .& (data[!, "DEC"].>-3.0)),:]
    G15 = data[((data[!, "RA"].<223.5) .& (data[!, "RA"].>211.5) .& (data[!, "DEC"].<3.0) .& (data[!, "DEC"].>-2.0)),:]
    G23 = data[((data[!, "RA"].<351.9) .& (data[!, "RA"].>338.1) .& (data[!, "DEC"].<-30.0) .& (data[!, "DEC"].> -35.0)),:]

    println(findHighestDensityArea(G02, 3))

    println("Done")
end


main()
