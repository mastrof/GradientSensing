##
using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames, Distributions, HypothesisTests

##
function phycosphere_kolmogorovsmirnov()
    datasets = collect_results(
        datadir("Poisson", "KolmogorovSmirnov");
        rinclude = [r"sensing"]
    ) |> unpack_dataframe
 
    dfgroups = [:T, :U, :Dc]
    gdf = groupby(datasets, dfgroups)

    for df in gdf
        C₀ = df.C₀[1]
        Dc = df.Dc[1]
        T = df.T[1]
        U = df.U[1]
        N = df.N[1]
        Δt = df.Δt[1]
        config = @strdict C₀ Dc T U N Δt
        produce_or_load(
            datadir("Poisson", "KolmogorovSmirnov"), config;
            prefix = "phycosphere",
            suffix = "jld2",
            tag = false,
            loadfile = false,
            force = true
        ) do config
            S = phycosphere_kolmogorovsmirnov(df)
            @strdict S
        end
    end
end

function phycosphere_kolmogorovsmirnov(df)
    f = jldopen(datadir("Poisson", "RC.jld2"), "r")
    R, Cₛ = f["R"], f["Cₛ"]
    close(f)
    iter = Iterators.product(R, Cₛ) |> collect
    S = zeros(typeof(1.0u"μm"), size(iter))
    println(size(S))
    for i in eachindex(iter)
        try
            R, Cₛ = iter[i]
            row = subset(df,
                :R => r -> r .== round(ustrip(R), sigdigits=3),
                :Cₛ => c -> c .== round(ustrip(Cₛ), sigdigits=3)
            ) |> first
            l = length(row.ksavg)
            r = [collect(s[end-l+1:end]) for s in row.r]
            mr = mean(r)
            S[i] = sensing_threshold(mr, row.ksavg, R)
        catch
            S[i] = NaN * 1u"μm"
        end
    end
    S
end

"""
    sensing_threshold(transect, sensing, R)
Finds the first position along the transect where sensing
(as defined by the average of the one-tailed KS test over independent samples)
occurs with a probability greater than 0.99.

The transect is in reverse order, i.e. it starts away from the source and
moves towards it. If there is no position at which `sensing > 0.99`,
then the function returns the radius of the source (`R`) itself,
meaning that there is no distance at which sensing can occur consistently.
"""
function sensing_threshold(grid, sensing, R)
    idx = findfirst(sensing .> 0.99)
    isnothing(idx) ? R : grid[idx]
end

##
phycosphere_kolmogorovsmirnov()
