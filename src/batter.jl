function parse_batter(batter_output)
    df = CSV.read(batter_output, DataFrame; header = ["name", "start", "end", "id", "score", "strand", "fppk"])
    return [(r.start, r.strand[1]) for r in eachrow(df)]
end

function batter(input; kwargs...)
    if haskey(kwargs, :batter_output)
        return batter(input, output)
    else
        output = tempname()
        results = batter(input, output)
        rm(output)
        return results
    end
end

function batter(genome::GenomicAnnotations.Record, output; kwargs...)
    tmpfasta = tempname()
    open(FASTA.Writer, tmpfasta) do fna
        write(fna, FASTA.Record("genome", genome.sequence))
    end
    batter(tmpfasta, output; kwargs...)
end


function batter(fasta::AbstractString, output; kwargs...)
    args = join([string("--", k, " ", v) for (k, v) in kwargs], " ")
    if TranscriptionUnits.cuda_available && !isempty(args)
        Conda.runconda(`run -n $batter_env $batter_tpe --fasta $fasta --output $output $args`, batter_env)
    elseif TranscriptionUnits.cuda_available && isempty(args)
        Conda.runconda(`run -n $batter_env $batter_tpe --fasta $fasta --output $output`, batter_env)
    elseif !TranscriptionUnits.cuda_available && !isempty(args) && !haskey(kwargs, "device")
        Conda.runconda(`run -n $batter_env $batter_tpe --fasta $fasta --output $output --device cpu $args`, batter_env)
    elseif !TranscriptionUnits.cuda_available && !isempty(args)
        Conda.runconda(`run -n $batter_env $batter_tpe --fasta $fasta --output $output $args`, batter_env)
    else
        Conda.runconda(`run -n $batter_env $batter_tpe --fasta $fasta --output $output --device cpu`, batter_env)
    end
    return parse_batter(output)
end
