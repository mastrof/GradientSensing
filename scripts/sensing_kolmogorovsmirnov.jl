using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames, Distributions, HypothesisTests

function sensing_kolmogorovsmirnov()
    f = jldopen(datadir("Poisson", "RC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]
    sensing_from_waitingtimes(R, Cₛ)
end

function sensing_from_waitingtimes(R::AbstractVector, Cₛ::AbstractVector)
    for i in eachindex(R)
        datasets = get_datasets_waitingtimes(R[i])
        foreach(sensing_from_waitingtimes, eachrow(datasets))
    end
end

"""
    get_datasets_waitingtimes(R)
Collects all the waiting time datasets for source radius `R`.
"""
function get_datasets_waitingtimes(R)
    Rsig3 = round(typeof(1.0u"μm"), R, sigdigits=3) |> ustrip
    collect_results(
        datadir("Poisson");
        rinclude = [Regex("waitingtimes.*R=$(Rsig3)")]
    )
end

"""
    sensing_from_waitingtimes(df)
Performs a Kolmogorov-Smirnov sensing test over the data stored in
the dataframe `df`. See `kstest_transect`.
Outputs (transect midpoints, KS sensing statistics, fraction of KS successes)
are saved to file (with suffix "sensing").
"""
function sensing_from_waitingtimes(df)
    params = parse_savename(df.path)[2]
    produce_or_load(
        datadir("Poisson", "KolmogorovSmirnov"), params;
        prefix="sensing", suffix="jld2", tag=false, loadfile=false
    ) do params
        r = df.r
        waitingtimes = df.waitingtimes
        R = params["R"]u"μm"
        Cₛ = params["Cₛ"]u"μM"
        Δt = params["Δt"]u"ms"
        ks = kstest_transect(waitingtimes, r, Δt, R, Cₛ, C₀)
        # evaluate mean in each interval
        ksavg = mean(ks)
        @strdict r ks ksavg
    end
end

# aliases for brevity
VVQ = AbstractVector{<:AbstractVector{<:Quantity}}
VVVQ = AbstractVector{<:AbstractVector{<:AbstractVector{<:Quantity}}}
function kstest_transect(waitingtimes::VVVQ, r, Δt, R, Cₛ, C₀)
    [kstest_transect(w,r,Δt,R,Cₛ,C₀) for w in waitingtimes]
end
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

sensing_kolmogorovsmirnov()
