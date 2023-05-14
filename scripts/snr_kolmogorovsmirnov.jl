using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames

function snr_at_ksphycosphere()
    datasets = collect_results(
        datadir("Poisson", "KolmogorovSmirnov"),
        rinclude = [r"phycosphere"]
    )
    for df in eachrow(datasets)
        snr_at_phycosphere(df)
    end
end

function snr_at_phycosphere(df)
    f = jldopen(datadir("Poisson", "RC.jld2"))
    @unpack R, Cₛ = f
    params = parse_savename(df.path)[2]
    params["Π"] = 1
    produce_or_load(
        datadir("Poisson", "KolmogorovSmirnov"), params;
        prefix = "SNR", suffix = "jld2",
        tag = false, loadfile = false, force = true
    ) do params
        S = df.S # size(S) == (length(R), length(Cₛ))
        for j in axes(S,2), i in axes(S,1)
            if S[i,j] == R[i]
                S[i,j] = NaN * unit(S[i,j])
            end
        end
        T = params["T"]u"ms"
        U = params["U"]u"μm/s"
        Π = params["Π"]u"1"
        SNR = @. snr(S, R, Cₛ', C₀, T, U, Π)
        @strdict SNR
    end
end

snr_at_ksphycosphere()
