using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames

function ic_hein()
    f = jldopen(datadir("Hein", "paramspaceRL.jld2"))
    R, L = f["R"], f["L"]

    # run length used to evaluate IC
    λ = 46.5/2#u"μm"

    datasets = collect_results(datadir("Hein"), rinclude=[r"phycosphere"])
    for df in eachrow(datasets)
        params = parse_savename(df.path)[2]
        # add λ to params
        params["λ"] = λ
        S = df.S

        produce_or_load(
            datadir("Hein"), params;
            prefix="IC", suffix="jld2",
            tag=false
        ) do params
            # reassign units to parameters
            λ = params["λ"]u"μm"
            # evaluate chemotactic index for each (R,L) pair
            ic = @. IC(λ,R,S)
            @strdict ic
        end
    end
end

ic_hein()