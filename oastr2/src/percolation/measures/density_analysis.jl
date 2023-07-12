#TO DO: in use but must be updated to recieve MST inputs.

using Random, StatsBase

L = 100
n_sims = 100

rho_min = 0.1
rho_max = 0.9
rho_step = 0.01
rho_range = rho_min:rho_step:rho_max

function percolation_cluster(rho)
    grid = rand(Float64, L, L)
    cluster = zeros(Int, L, L)
    for i in 1:L, j in 1:L
        if grid[i, j] < rho
            cluster[i, j] = 1
        end
    end
    labels = label_components(cluster)
    return maximum(labels)
end

function newman_ziff(L, rho)
    n_clusters = zeros(Int, n_sims)
    for i in 1:n_sims
        n_clusters[i] = percolation_cluster(rho)
    end
    p_inf = count(x -> x == L^2, n_clusters) / n_sims
    return p_inf
end

for rho in rho_range
    p_inf = newman_ziff(L, rho)
    println("$rho $p_inf")
end
