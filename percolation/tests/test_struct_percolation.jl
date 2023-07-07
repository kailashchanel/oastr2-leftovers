#test of struct features from Newmann-Ziff Algorithm Paper. See paper draft for explaination. Not in use.

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


using Test

G1 = QuickFind(10);
G2 = QuickUnion(10, false, false);
G3 = QuickUnion(10, true, false);
G4 = QuickUnion(10, false, true);
G5 = QuickUnion(10, true, true);

function test(G)

    connect!(G, 1, 6)
    connect!(G, 6, 7)
    connect!(G, 7, 2)
    connect!(G, 2, 3)
    connect!(G, 3, 8)
    
    connect!(G, 9, 4)
    connect!(G, 4, 5)
    connect!(G, 5, 10)

    # `pts1` should form one set of connected components and `pts2` another, respectively.
    pts1 = [1, 2, 3, 6, 7, 8]
    pts2 = [4, 5, 9, 10]
    for p in pts1
        for q in pts1
            @test isconnected(G, p, q)
            @test isconnected(G, q, p)
        end
        for q in pts2
            @test !isconnected(G, p, q)
            @test !isconnected(G, q, p)
        end
    end
end

@testset "Test" begin
    @testset "Test QuickFind" begin test(G1) end
    @testset "Test QuickUnion without weighting and path compression" begin test(G2) end
    @testset "Test QuickUnion with weighting and without path compression" begin test(G3) end
    @testset "Test QuickUnion without weighting and with path compression" begin test(G4) end
    @testset "Test QuickUnion with weighting and path compression" begin test(G5) end
end
