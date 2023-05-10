using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames, Distributions, HypothesisTests

"""
    kstestright(waitingtimes, c, Δt)
Performs a one-tailed Kolmogorov-Smirnov test between sampled `waitingtimes`
and a reference distribution of the absorption waiting times that
would be observed at a constant concentration `c`.
The reference distribution is a truncated exponential with lower bound `Δt`
and average `waitingtime(c)`.

The right-tailed Kolmogorov-Smirnov tests the following null hypothesis:
**the true CDF of the sample is not greater than the CDF of the reference distribution**.
The hypothesis is tested by first evaluating the KS statistics with the
`ExactOneSampleKSTest` function, and then evaluating the pvalue of the right tail
of this statistics (via `pvalue(KStest, tails=:right)`).
In this scenario, a low p-value implies a high likelihood that the true CDF of
the sample is larger than (i.e. it lies above) the CDF of the reference distribution.

Therefore, if p<0.05 the sample shows significantly shorter waiting times than
the reference distribution.
If this condition is satisfied, the function returns `1.0`
(i.e. the sensor has observed a statistically significant decrease in waiting times),
otherwise (p≥0.05) the function returns `0.0`.
"""
function kstestright(waitingtimes, c, Δt)
    isempty(waitingtimes) && return Float64(false)
    # need to convert all values to same units and then ustrip
    # because HypothesisTests functions don't work with units
    τ = waitingtime(c) |> u"ms" |> ustrip
    tmin = Δt |> u"ms" |> ustrip
    # expected pdf of waiting times at concentration c
    reference_pdf = Truncated(Exponential(τ), tmin, Inf)
    # compare sampled waiting times to reference pdf via Kolmogorov-Smirnov test
    KStest = ExactOneSampleKSTest(ustrip.(waitingtimes), reference_pdf)
    # if the right-tailed p < 0.05, the empirical CDF is larger than the reference
    return Float64(pvalue(KStest, tail=:right) < 0.05)
end

# aliases for brevity
VVQ = AbstractVector{<:AbstractVector{<:Quantity}}
VVVQ = AbstractVector{<:AbstractVector{<:AbstractVector{<:Quantity}}}
"""
    kstest_transect(waitingtimes, r, Δt, R, Cₛ, C₀)
Performs a Kolmogorov-Smirnov test (see `kstestright`) for the `waitingtimes`
at each point in the transect described by `r`.
`Δt` is the temporal resolution which the sensor has sampled the absorption events.
`R`, `Cₛ` and `C₀` are used to evaluate the real concentration field at the midpoint
of each interval in the transect, which defines the reference distribution for
the absorption waiting times in the KS test.
"""
function kstest_transect(waitingtimes::VVQ, r, Δt, R, Cₛ, C₀)
    c = @. C(r, R, Cₛ, C₀)
    [kstestright(waitingtimes[i], c[i], Δt) for i in eachindex(c)]
end
function kstest_transect(waitingtimes::VVVQ, r, Δt, R, Cₛ, C₀)
    [kstest_transect(w,r,Δt,R,Cₛ,C₀) for w in waitingtimes]
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

function phycosphere_ks()
    f = jldopen(datadir("PoissonSampling", "paramspaceRC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]

    # loop over R to only read a few files at once
    # otherwise memory would get overloaded
    for i in eachindex(R)
        Rᵢ = round(ustrip(R[i]), sigdigits=3)
        datasets = collect_results(datadir("PoissonSampling");
            rinclude=[Regex("waitingtimes.*R=$(Rᵢ)_")]
        )
        for df in eachrow(datasets)
            params = parse_savename(df.path)[2]
            produce_or_load(
                datadir("PoissonSampling", "KolmogorovSmirnov"), params;
                prefix="sensing", suffix="jld2", tag=false, loadfile=false
            ) do params
                r = df.r
                waitingtimes = df.waitingtimes
                # C₀ = params["C₀"]u"nM"
                local R = params["R"]u"μm"
                local Cₛ = params["Cₛ"]u"μM"
                Δt = params["Δt"]u"ms"
                ks = kstest_transect(waitingtimes, r, Δt, R, Cₛ, C₀)
                # evaluate mean in each interval
                ksavg = mean(ks)
                @strdict r ks ksavg
            end
        end
    end

    # set parameters
    T = 100u"ms"
    Δt = 1e-6u"ms"
    U = 46.5u"μm/s"
    N = 100
    params = @strdict C₀ T U Δt N
    params_ustrip = Dict(keys(params) .=> ustrip.(values(params)))
    produce_or_load(
        datadir("PoissonSampling", "KolmogorovSmirnov"), params_ustrip;
        prefix="phycosphere", suffix="jld2", tag=false, force=true
    ) do params_ustrip
        iter = Iterators.product(R,Cₛ) |> collect
        S = zeros(typeof(1.0u"μm"), size(iter))
        for i in eachindex(iter)
            R, Cₛ = iter[i]
            p = Dict("R" => ustrip(R), "Cₛ" => ustrip(Cₛ), params_ustrip...)
            fname = datadir("PoissonSampling", "KolmogorovSmirnov",
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


phycosphere_ks()