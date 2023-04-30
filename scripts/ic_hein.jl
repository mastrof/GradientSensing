using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames

"""
Evaluates the black-hole chemotactic index for the Hein phycosphere.

Reads the value of phycosphere radii (`S`) evaluated for each `(R,Cₛ)` pair
by `scripts/phycosphere_hein.jl` and evaluates the IC using the formulas
defined in `src/chemotactic_index.jl`.

Returns a jld2 file with the IC value for each `(R,Cₛ)` pair.
"""
function ic_hein()
    f = jldopen(datadir("Hein", "paramspaceRC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]

    # run length used to evaluate IC
    λ = 46.5u"μm"/2

    datasets = collect_results(datadir("Hein"), rinclude=[r"phycosphere"])
    for df in eachrow(datasets)
        params = parse_savename(df.path)[2]
        # add λ to params
        params["λ"] = ustrip(λ)
        S = df.S

        produce_or_load(
            datadir("Hein"), params;
            prefix="IC", suffix="jld2",
            tag=false
        ) do params
            # evaluate chemotactic index for each (R,L) pair
            ic = @. IC(λ,R,S)
            @strdict ic
        end
    end
end

ic_hein()