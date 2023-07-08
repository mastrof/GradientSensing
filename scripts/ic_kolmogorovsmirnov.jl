using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames

function ic_ks()
    f = jldopen(datadir("Poisson", "RC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]

    # run length used to evaluate IC
    λ = 30u"μm"

    datasets = collect_results(
        datadir("Poisson", "KolmogorovSmirnov"),
        rinclude=[r"phycosphere"]
    )

    for df in eachrow(datasets)
        params = parse_savename(df.path)[2]
        # add λ to params
        params["λ"] = ustrip(λ)
        S = df.S

        produce_or_load(
            datadir("Poisson", "KolmogorovSmirnov"), params;
            prefix="IC", suffix="jld2", tag=false
        ) do params
            ic = @. IC(λ, R, S)
            @strdict ic
        end
    end
end

ic_ks()
