using TranscriptionUnits
using GenomicAnnotations
using Test

@testset "TranscriptionUnits" begin
    gbk = readgbk("A0901.gbk")[1]
    boundaries = transcriptionunits(gbk, TranscriptionUnits.batter(gbk))
    @test !isempty(boundaries)
end