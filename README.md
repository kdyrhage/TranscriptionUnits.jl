# TranscriptionUnits.jl
Predictor for transcriptional units.

## Installation
While running Julia interactively, press `]` to enter the package manager and run:
```julia
add https://github.com/kdyrhage/TranscriptionUnits.jl
```

## Usage
The package provides the function `transcriptionunits(genome)`. `genome` can be the path to a file containing annotations in the GenBank format, or a genomic `Record` from [`GenomicAnnotations`](https://github.com/BioJulia/GenomicAnnotations.jl). The function also accepts the following optional keyword arguments:
- `distthreshold`: maximum number of nucleotides allowed between genes in a unit. (Default: 7)
- `terminators`: an array where each element is a `Tuple{Int, Char}` containing start positions and strand ('+'/'-') for predicted transcriptional terminators. This is only useful for large `distthreshold`s. The package incorporates a predictor, [BATTER](https://github.com/lulab/BATTER), which can be called using `batter(chr)`. For best performance BATTER should be run on a system with a functional CUDA installation, but if none can be detected it will be executed using the CPU instead.

## Example
```julia
using TranscriptionUnits
TUs = transcriptionunits("test/A0901.gbk")
open("output.txt", "w") do out
    println(out, first(TUs))
end
```

```julia
using GenomicAnnotations
chr = readgbk("test/A0901.gbk")[1]
TUs = transcriptionunits(chr, distthreshold = 100, terminators = batter(chr))

genes = @genes(chr, CDS)

for TU in TUs
    println(join(genes[TU].locus_tag, ","))
end
```