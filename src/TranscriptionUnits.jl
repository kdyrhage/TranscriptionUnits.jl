module TranscriptionUnits

using DataFrames
using Scratch
using CSV
using GenomicAnnotations
using FASTX
using Conda

batter_root = ""
const batter_env = :batter_env
batter_tpe = ""
cuda_available = false

const defaults = Dict(:overlap => 0)

function get_batter()
    global batter_root = @get_scratch!("batter")
    if isempty(readdir(batter_root))
        @eval begin
            using Git
            run(`$(git()) clone https://github.com/lulab/BATTER.git $batter_root`)
        end
    end
    prerequisites = ["pyfaidx", "transformers", "pytorch", "pandas", "numpy"]
    if length(intersect(prerequisites, Conda._installed_packages(batter_env))) != length(prerequisites)
        condachannels = ["pytorch", "bioconda", "conda-forge", "defaults"]
        if Conda.channels(batter_env) != condachannels
            foreach(c -> Conda.add_channel(c, batter_env), reverse(condachannels))
        end
        Conda.add(["pyfaidx=0.7.1", "transformers=4.18", "pytorch=1.8", "pandas", "numpy=1", "mkl=2020"], batter_env)
    end
    return batter_root
end

function __init__()
    batter_root = get_batter()
    python = joinpath(Conda.python_dir(:batter_env), "python")
    out = IOBuffer()
    run(pipeline(ignorestatus(Conda._set_conda_env(`$(python) -c "import torch;print(torch.cuda.is_available())"`, :batter_env)), stdout = out, stderr = devnull))
    global cuda_available = String(take!(out)) == "True\n"
    global batter_tpe = joinpath(batter_root, "scripts", "batter-tpe")
end

include("batter.jl")
include("predictor.jl")

export transcriptionunits

end # module TranscriptionalUnits
