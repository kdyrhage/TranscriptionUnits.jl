struct Directone
    first::Int
    last::Int
end

genedist(g::GenomicAnnotations.Record) = genedist(g.genes)
genedist(g::AbstractVector{Gene}) = [locus(g[i]).start - locus(g[i - 1]).stop - 1 for i in 2:size(g, 1)]
genedist(g::AbstractVector{Gene}, directone::Directone) = genedist(g, directone.first, directone.last)
function genedist(g::AbstractVector{Gene}, firstgene::Int, lastgene::Int)
    firstgene == lastgene && return Int[]
    Int[locus(g[i]).start - locus(g[i - 1]).stop - 1 for i in (firstgene + 1):lastgene]
end

directones(genome::GenomicAnnotations.Record) = directones(genome.genes)
function directones(g::AbstractVector{Gene})
    boundaries = Directone[]
    strand = iscomplement(g[1])
    directonestart = 1
    for i in 2:size(g, 1)
        if iscomplement(g[i]) != strand
            strand = iscomplement(g[i])
            push!(boundaries, Directone(directonestart, i - 1))
            directonestart = i
        end
    end
    push!(boundaries, Directone(directonestart, size(g, 1)))
    boundaries
end

function tu_boundaries(directones::AbstractVector{Directone}, dists::AbstractVector{Vector{Int}}, threshold = 55)
    boundaries = []
    for d in 1:length(directones)
        first = directones[d].first
        last = directones[d].last
        append!(boundaries, [first, last])
        @assert (last - first) == length(dists[d])
        for i in 1:length(dists[d])
            if dists[d][i] > threshold
                append!(boundaries, [first + i - 1, first + i])
            end
        end
    end
    sort!(boundaries)
    [(boundaries[i], boundaries[i+1]) for i in 1:2:length(boundaries)]
end

function isintergenic(g::AbstractVector{Gene}, pos)
    for i in eachindex(g)
        locus(g[i]).start <= pos <= locus(g[i]).stop && feature(g[i]) in [:gene, :CDS, :tRNA, :rRNA] && return false
    end
    true
end

function isintergenic(g::AbstractVector{Gene}, pos, o)
    for i in eachindex(g)
        overlap = o * (locus(g[i]).stop - locus(g[i]).start)
        (locus(g[i]).start + overlap) <= pos <= (locus(g[i]).stop - overlap) && feature(g[i]) in [:gene, :CDS, :tRNA, :rRNA] && return false
    end
    true
end

function within(boundaries::Vector{Tuple{Int, Int}}, pos)
    for i in eachindex(boundaries)
        minimum(boundaries[i]) <= pos <= maximum(boundaries[i]) && return i
    end
    0
end

function containsgenes(g::AbstractVector{Gene}, boundaries::Tuple{Int, Int})
    genecount = 0
    for i in eachindex(g)
        if locus(g[i]).start >= first(boundaries)
            locus(g[i]).stop <= last(boundaries) ? (genecount += 1) : return genecount
        end
    end
    return genecount
end

function addbreak!(g::AbstractVector{Gene}, boundaries::Vector{Tuple{Int, Int}}, t, i)
    boundary = boundaries[i]
    deleteat!(boundaries, i)
    firsthalf = (first(boundary), t - 1)
    lasthalf = (t, last(boundary))
    containsgenes(g, firsthalf) >= 2 && push!(boundaries, (first(boundary), t - 1))
    containsgenes(g, lasthalf) >= 2 && push!(boundaries, (t, last(boundary)))
end

function addterminators!(g::AbstractVector{Gene}, boundaries::Vector{Tuple{Int64, Int64}}, terminators::Vector{Int}, overlap = 0)
    for t in terminators
        i = within(boundaries, t)
        i > 0 && isintergenic(g, t, overlap) && addbreak!(g, boundaries, t, i)
    end
    sort!(boundaries)
    nothing
end

function issamestrand(genes::AbstractVector{Gene}, t::Tuple{Int64, Char})
    termindatorstrand = t[2] == '+'
    for i in 2:length(genes)
        if locus(genes[i - 1]).stop <= t[1] <= locus(genes[i]).start
            iscomplement(genes[i - 1]) == iscomplement(genes[i]) == terminatorstrand && return true
        end
    end
    return false
end

function addterminators!(g::AbstractVector{Gene}, boundaries::Vector{Tuple{Int64, Int64}}, terminators::Vector{Tuple{Int64, Char}}, overlap = 0)
    for t in terminators
        i = within(boundaries, t[1])
        i > 0 && isintergenic(g, t[1], overlap) && issamestrand(g, t) && addbreak!(g, boundaries, t[1], i)
    end
    sort!(boundaries)
    nothing
end


transcriptionunits(g, d = 7, t = Tuple{Int, Char}[]; overlap = 0, batter_output = tempname()) =
    transcriptionunit(g, d, t; kwargs...)
transcriptionunits(g::GenomicAnnotations.Record, d::Function, t; kwargs...) =
    transcriptionunits(g, d(genedist(g.genes)), t; kwargs...)
transcriptionunits(genomes::AbstractVector{GenomicAnnotations.Record}, d, t; kwargs...) =
    map(g -> transcriptionunits(g, d, t; kwargs...), genomes)
transcriptionunits(g::GenomicAnnotations.Record, d::Int; kwargs...) = transcriptionunits(g, d, batter(g; kwargs...); kwargs...)
function transcriptionunits(genome::GenomicAnnotations.Record, distthreshold, terminators; kwargs...)
    d = directones(genome)
    dists = [genedist(genome.genes, directone) for directone in d]
    boundaries = tu_boundaries(d, dists, distthreshold)
    addterminators!(genome.genes, boundaries, terminators, get(kwargs, :overlap, 0))
    return boundaries
end