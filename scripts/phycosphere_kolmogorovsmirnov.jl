using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames, Distributions, HypothesisTests

function phycosphere_kolmogorovsmirnov()
    f = jldopen(datadir("newPoisson", "RC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]
    R = R[1:end-1]

    # set parameters
    T = 100u"ms"
    Δt = 1e-4u"ms"
    U = 46.5u"μm/s"
    N = 250
    params = @strdict C₀ T U Δt N
    params_ustrip = Dict(keys(params) .=> ustrip.(values(params)))
    produce_or_load(
        datadir("newPoisson", "KolmogorovSmirnov"), params_ustrip;
        prefix="phycosphere", suffix="jld2", tag=false
    ) do params_ustrip
        iter = Iterators.product(R,Cₛ) |> collect
        S = zeros(typeof(1.0u"μm"), size(iter))
        for i in eachindex(iter)
            R, Cₛ = iter[i]
            p = Dict("R" => ustrip(R), "Cₛ" => ustrip(Cₛ), params_ustrip...)
            fname = datadir("newPoisson", "KolmogorovSmirnov",
                savename("sensing", p, "jld2")
            )
            @unpack r, ks, ksavg = load(fname)
            #S[i] = mean(sensing_threshold.(r, ks, R))
            mr = mean(collect.(r))
            S[i] = sensing_threshold(mr, ksavg, R)
        end
        @strdict S
    end
end

"""
    sensing_threshold(transect, sensing, R)
Finds the first position along the transect where sensing
(as defined by the average of the one-tailed KS test over independent samples)
occurs with a probability greater than chance (>0.5).

The transect is in reverse order, i.e. it starts away from the source and
moves towards it. If there is no position at which `sensing > 0.5`,
then the function returns the radius of the source (`R`) itself,
meaning that there is no distance at which sensing can occur consistently.
"""
function sensing_threshold(grid, sensing, R)
    idx = findfirst(sensing .> 0.5)
    isnothing(idx) ? R : grid[idx]
end

phycosphere_kolmogorovsmirnov()
