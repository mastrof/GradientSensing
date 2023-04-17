using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames

function ic_ks()
    f = jldopen(datadir("PoissonSampling", "paramspaceRL.jld2"))
    R, L = f["R"], f["L"]

    # run length used to evaluate IC
    λ = 46.5u"μm"/2

    datasets = collect_results(
        datadir("PoissonSampling", "KolmogorovSmirnov"),
        rinclude=[r"phycosphere"]
    )

    for df in eachrow(datasets)
        params = parse_savename(df.path)[2]
        # add λ to params
        params["λ"] = ustrip(λ)
        S = df.S

        produce_or_load(
            datadir("PoissonSampling", "KolmogorovSmirnov"), params;
            prefix="IC", suffix="jld2", tag=false
        ) do params
            ic = @. IC(λ, R, S)
            @strdict ic
        end
    end
end

ic_ks()