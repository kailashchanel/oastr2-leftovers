#percolation test on a latice. Not in use.

using Random

function percolate_lattice(n::Int, p::Float64)
    # initialize the lattice with all bonds intact
    lattice = ones(Int8, n, n)

    # generate a random number generator
    rng = Random.default_rng()

    # randomly remove bonds with probability p
    for i in 1:n-1, j in 1:n-1
        if rng.uniform() < p
            # remove the right bond
            lattice[i, j] &= 0
            lattice[i+1, j] &= 0
        end
        if rng.uniform() < p
            # remove the bottom bond
            lattice[i, j] &= 0
            lattice[i, j+1] &= 0
        end
    end

    return lattice
end

function percolate_cluster(lattice::Array{Int8, 2}, i::Int, j::Int)
    # check if the bond at (i,j) is intact
    if lattice[i, j] == 0
        return
    end

    # mark the bond as visited
    lattice[i, j] = 0

    # recursively visit neighboring bonds
    if i > 1
        percolate_cluster(lattice, i-1, j)
    end
    if i < size(lattice, 1)
        percolate_cluster(lattice, i+1, j)
    end
    if j > 1
        percolate_cluster(lattice, i, j-1)
    end
    if j < size(lattice, 2)
        percolate_cluster(lattice, i, j+1)
    end
end

function percolate_probability(n::Int, p::Float64, trials::Int)
    percolation_probabilities = zeros(Float64, trials)

    for trial in 1:trials
        # generate a percolated lattice
        lattice = percolate_lattice(n, p)

        # check if the top and bottom of the lattice are connected
        percolates = false
        for j in 1:n
            if lattice[1, j] == 1
                percolate_cluster(lattice, 1, j)
            end
            if lattice[n, j] == 1
                percolate_cluster(lattice, n, j)
            end
            if lattice[1, j] == 0 && lattice[n, j] == 0
                percolates = true
                break
            end
        end

        # record the percolation probability
        percolation_probabilities[trial] = percolates
    end

    return mean(percolation_probabilities)
end
