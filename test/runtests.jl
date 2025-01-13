using TranscriptionUnits
using GenomicAnnotations
using Conda
using Test

@testset "TranscriptionUnits" begin
    chr = readgbk("A0901.gbk")[1]
    boundaries = transcriptionunits(chr)
    @test !isempty(boundaries)
end

#@testset "BATTER" begin
#    Conda.runconda(`run -n $(TranscriptionUnits.batter_env) $(TranscriptionUnits.batter_tpe) -h`, TranscriptionUnits.batter_env)
#end
