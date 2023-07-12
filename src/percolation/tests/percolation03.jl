#Not in use. This is an alternate method using a direct implemenation of the Newmann-Ziff Algorithm

import Pkg
Pkg.add(["JLD2, Plots"])

using JLD2, Plots

struct QuickFind
    id :: Vector{Int}
end

function QuickFind(n::Int)
    id = collect(1:n)
    return QuickFind(id)
end

function find(qf::QuickFind, p::Int)
    id = qf.id
    return id[p]
end

function isconnected(qf::QuickFind, p::Int, q::Int)
    i = find(qf, p)
    j = find(qf, q)
    return isequal(i, j)
end

function connect!(qf::QuickFind, p::Int, q::Int)
    id = qf.id
    pid = id[p]
    qid = id[q]
    n = length(id)
    for i=1:n
        if isequal(getindex(id, i), pid)
            setindex!(id, qid, i)
        end
    end
    return
end

"""
    QuickUnion
"""
struct QuickUnion
    id :: Vector{Int}
    sz :: Vector{Int}
    weighted :: Bool
    compress :: Bool
end

function QuickUnion(n::Int, weighted=true, compress=true)
    id = collect(1:n)
    sz = zeros(Int, n)
    return QuickUnion(id, sz, weighted, compress)
end

function find(qu::QuickUnion, p::Int)
    id = qu.id
    while p != id[p]
        if qu.compress
            id[p] = id[id[p]]
        end
        p = id[p]
    end
    return p
end

function isconnected(qu::QuickUnion, p::Int, q::Int)
    i = find(qu, p)
    j = find(qu, q)
    return isequal(i, j)
end

"""
    connect(G, p, q)
Connect p and q in G.
If weighting is used, connect such a way that the depth of the tree is minimized.
"""
function connect!(qu::QuickUnion, p::Int, q::Int)
    id = qu.id
    i = find(qu, p)
    j = find(qu, q)
    if qu.weighted
        sz = qu.sz
        if (sz[i] < sz[j])
            id[i] = j
            sz[j] += sz[i]
        else
            id[j] = i
            sz[i] += sz[j]
        end
    else
        id[i] = j
    end
end

export QuickFind, QuickUnion, connect!, isconnected, find

struct Percolation
    grid :: Matrix{Bool}
    li :: LinearIndices
    ny :: Int
    nx :: Int
    wuf1 :: QuickUnion
    wuf2 :: QuickUnion
    node1 :: Int
    node2 :: Int
end

function Percolation(ny, nx)
    grid = zeros(Bool, ny, nx)
    li = LinearIndices(grid)
    wuf1 = QuickUnion(nx*ny + 2)
    wuf2 = QuickUnion(nx*ny + 1)
    return Percolation(grid, li, ny, nx, wuf1, wuf2, nx*ny + 1, nx*ny + 2)
end

function connect!(p::Percolation, i, j)
    connect!(p.wuf1, i, j)
    connect!(p.wuf2, i, j)
end

function isopen(p::Percolation, i, j)
    return p.grid[i, j]
end

function number_of_open(p::Percolation)
    return sum(p.grid)
end

function isfull(p::Percolation, i, j)
    !isopen(p, i, j) && return false
    return isconnected(p.wuf2, p.node1, p.li[i, j])
end

function percolates(p::Percolation)
    return isconnected(p.wuf1, p.node1, p.node2)
end

function open!(p::Percolation, i, j)
    isopen(p, i, j) && return
    p.grid[i, j] = true

    i == 1       && connect!(p, p.node1, p.li[i, j])
    i == p.ny    && connect!(p.wuf1, p.node2, p.li[i, j])

    i > 1        && isopen(p, i-1, j) && connect!(p, p.li[i, j], p.li[i-1, j])
    i < p.ny - 1 && isopen(p, i+1, j) && connect!(p, p.li[i, j], p.li[i+1, j])
    j > 1        && isopen(p, i, j-1) && connect!(p, p.li[i, j], p.li[i, j-1])
    j < p.nx - 1 && isopen(p, i, j+1) && connect!(p, p.li[i, j], p.li[i, j+1])
end

using Statistics

struct Result{T}
    samples :: Vector{T}
    mean :: Float64
    stddev :: Float64
    confidence_lo :: Float64
    confidence_hi :: Float64
end

const CONFIDENCE_95 = 1.96

function Result(samples)
    n = length(samples)
    m = mean(samples)
    s = std(samples)
    cl = m - CONFIDENCE_95 * s / sqrt(n)
    ch = m + CONFIDENCE_95 * s / sqrt(n)
    return Result(samples, m, s, cl, ch)
end

using Printf

function Base.show(io::IO, r::Result)
    @printf("% 25s =  %d\n", "samples", length(r.samples))
    @printf("% 25s =  %0.3f\n", "mean", r.mean)
    @printf("% 25s =  %0.3f\n", "stddev", r.stddev)
    @printf("% 25s = [%0.3f, %0.3f]\n", "95% confidence interval", r.confidence_lo, r.confidence_hi)
end

function run_one(ny, nx)
    p = Percolation(ny, nx)
    while !percolates(p)
        open!(p, rand(1:ny), rand(1:nx))
    end
    return p
end

function run(n, trials)
    thresholds = zeros(trials)
    for i=1:trials
        p = run_one(n, n)
        thresholds[i] = number_of_open(p) / n^2
    end
    return Result(thresholds)
end

struct GalaxyDistance
    idx1::Int32
    idx2::Int32
    distance::Float32
end

sorted_distance = load_object("sorted_distance_list.jld2")

side = Int(-.5 + sqrt(0.25+2*length(sorted_distance)))

function galaxy_run(sorted_distance,n)
    wuf1 = QuickUnion(n)
    clusters = Int32[]
    for gd in sorted_distance
        connect!(wuf1, Int64(gd.idx1), Int64(gd.idx2))
        numclusters = length(Set(wuf1.id))
        #println(gd.idx1," ",gd.idx2," ",gd.distance," ",numclusters)
        push!(clusters, Int32(numclusters))
    end
    return clusters
end

clusters = galaxy_run(sorted_distance,side)
Collapse
