using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames

function ic_ks()
    f = jldopen(datadir("newPoisson", "RC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]
    R = R[1:end-1]

    # run length used to evaluate IC
    λ = 46.5u"μm"/2

    datasets = collect_results(
        datadir("newPoisson", "KolmogorovSmirnov"),
        rinclude=[r"phycosphere"]
    )

    for df in eachrow(datasets)
        params = parse_savename(df.path)[2]
        # add λ to params
        params["λ"] = ustrip(λ)
        S = df.S

        produce_or_load(
            datadir("newPoisson", "KolmogorovSmirnov"), params;
            prefix="IC", suffix="jld2", tag=false, force=true
        ) do params
            ic = @. IC(λ, R, S)
            @strdict ic
        end
    end
end

ic_ks()
