using Distributed
@everywhere using DrWatson
@everywhere @quickactivate :GradientSensing
@everywhere begin
    using JLD2, DataFrames, Distributions, HypothesisTests
end


# aliases for brevity
@everywhere begin
    VVQ = AbstractVector{<:AbstractVector{<:Quantity}}
    VVVQ = AbstractVector{<:AbstractVector{<:AbstractVector{<:Quantity}}}
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

@everywhere function kstestright(waitingtimes, c, Δt, Dc)
    isempty(waitingtimes) && return Float64(false)
    # need to convert all values to same units and then ustrip
    # because HypothesisTests functions don't work with units
    τ = waitingtime(c, Dc) |> u"ms" |> ustrip
    tmin = Δt |> u"ms" |> ustrip
    # expected pdf of waiting times at concentration c
    reference_pdf = Truncated(Exponential(τ), tmin, Inf)
    # compare sampled waiting times to reference pdf via Kolmogorov-Smirnov test
    KStest = ExactOneSampleKSTest(ustrip.(waitingtimes), reference_pdf)
    # if the right-tailed p < 0.05, the empirical CDF is larger than the reference
    return Float64(pvalue(KStest, tail=:right) < 0.05)
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

@everywhere function kstest_transect(waitingtimes::VVQ, r, Δt, R, Cₛ, C₀, Dc)
    c = @. C(r, R, Cₛ, C₀)
    [kstestright(waitingtimes[i], c[i], Δt, Dc) for i in eachindex(waitingtimes)]
end
@everywhere function kstest_transect(waitingtimes::VVVQ, r, Δt, R, Cₛ, C₀, Dc)
    [kstest_transect(waitingtimes[i],r[i],Δt,R,Cₛ,C₀,Dc) for i in eachindex(waitingtimes)]
end


function get_filelist_waitingtimes()
    filenames = readdir(
        datadir("Poisson"); # on local
        #joinpath(ENV["SCRATCH"], "Poisson"); # on cluster
        join = true
    )
    filter!(s -> contains(s, "waitingtimes"), filenames)
    filenames
end

"""
    sensing_from_waitingtimes(df)
Performs a Kolmogorov-Smirnov sensing test over the data stored in
the dataframe `df`. See `kstest_transect`.
Outputs (transect midpoints, KS sensing statistics, fraction of KS successes)
are saved to file (with suffix "sensing").
"""

@everywhere function sensing_from_waitingtimes(fname)
    params = parse_savename(fname)[2]
    produce_or_load(
        datadir("Poisson", "KolmogorovSmirnov"), params; # on local
        #joinpath(ENV["SCRATCH"], "Poisson", "KolmogorovSmirnov"), params; # on cluster
        prefix="sensing", suffix="jld2", tag=false, loadfile=false
    ) do params
        f = jldopen(fname, "r")
        r = f["r"]
        waitingtimes = f["waitingtimes"]
        close(f)
        R = params["R"]u"μm"
        Cₛ = params["Cₛ"]u"μM"
        C₀ = params["C₀"]u"nM"
        Dc = params["Dc"]u"μm^2/s"
        Δt = params["Δt"]u"ms"
        ks = kstest_transect(waitingtimes, r, Δt, R, Cₛ, C₀, Dc)
        waitingtimes = nothing
        GC.gc()
        # evaluate mean in each interval
        l = minimum(length.(ks))
        #ksavg = mean([k[1:l] for k in ks])
        ksavg = mean([k[end-l+1:end] for k in ks])
        @strdict r ks ksavg
    end
    GC.gc()
    return nothing
end

"""
    produce_data(R::AbstractVector, Cₛ::AbstractVector)
Run sensing analysis on each dataset.
Collects datasets in chunks by values of R to avoid overloading memory.
"""
function produce_data()
    datasets = get_filelist_waitingtimes()
    pmap(sensing_from_waitingtimes, datasets)
    GC.gc()
end

##
produce_data()
