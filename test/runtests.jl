using TranscriptionUnits
using GenomicAnnotations
using Conda
using Test

@testset "TranscriptionUnits" begin
    chr = readgbk("A0901.gbk")[1]
    boundaries = transcriptionunits(chr)
    @test !isempty(boundaries)
end

@testset "BATTER" begin
    Conda.runconda(`run -n $batter_env $batter_tpe`, batter_envConda.runconda(`run -n $batter_env $batter_tpe -h`, batter_env))
end
