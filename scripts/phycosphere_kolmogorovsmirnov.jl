using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames, Distributions, HypothesisTests

function phycosphere_kolmogorovsmirnov()
    f = jldopen(datadir("Poisson", "RC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]

    # set parameters
    T = 100u"ms"
    Δt = 1e-4u"ms"
    U = 46.5u"μm/s"
    N = 100
    params = @strdict C₀ T U Δt N
    params_ustrip = Dict(keys(params) .=> ustrip.(values(params)))
    produce_or_load(
        datadir("Poisson", "KolmogorovSmirnov"), params_ustrip;
        prefix="phycosphere", suffix="jld2", tag=false
    ) do params_ustrip
        iter = Iterators.product(R,Cₛ) |> collect
        S = zeros(typeof(1.0u"μm"), size(iter))
        for i in eachindex(iter)
            R, Cₛ = iter[i]
            p = Dict("R" => ustrip(R), "Cₛ" => ustrip(Cₛ), params_ustrip...)
            fname = datadir("Poisson", "KolmogorovSmirnov",
                savename("sensing", p, "jld2")
            )
            # if the file is missing (e.g. not evaluated) assign NaNs
            @unpack r, ksavg = try
                load(fname)
            catch _;
                (r=NaN, ksavg=[NaN])
            end
            S[i] = isnan(ksavg[1]) ? (NaN)u"μm" : sensing_threshold(r, ksavg, R)
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
