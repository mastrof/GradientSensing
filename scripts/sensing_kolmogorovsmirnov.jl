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
    kstestright(waitingtimes, T)
Perform a one-tailed Kolmogorov-Smirnov test between the first and the
second half of the `waitingtimes` sampled within the interval of length `T`.

The right-tailed Kolmogorov-Smirnov tests the null hypothesis that
**the true CDF in the first half of the interval is not greater than
the CDF in the second half of the interval**.
In this setting, a low p-value indicates a high likelihood that the true
CDF of the sample in the first half of the interval is larger than the
CDF in the other half.

Therefore, if p<0.05, the first half shows significantly shorter
waiting times than the second half.
If this condition is satisfied, the function returns `1.0`
(i.e. the sensor has observed a statistically significant decrease
in waiting times), otherwise (p≥0.05) the function returns `0.0`.
"""
@everywhere function kstestright(waitingtimes, T)
    isempty(waitingtimes) && return Float64(false)
    # convert all values to same units and ustrip
    # because HypothesisTests does not support units
    dts = ustrip.(waitingtimes .|> u"ms")
    # split the sequence of waiting times around the midpoint T/2
    tsplit = ustrip(T/2 |> u"ms")
    isplit = findfirst(cumsum(dts) .> tsplit)
    dts_1 = @view dts[1:isplit-1]
    dts_2 = @view dts[isplit:end]
    KStest = ApproximateTwoSampleKSTest(dts_1, dts_2)
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

@everywhere function kstest_transect(waitingtimes::VVQ, r, R, Cₛ, C₀, T)
    c = @. C(r, R, Cₛ, C₀)
    [kstestright(waitingtimes[i], T) for i in eachindex(waitingtimes)]
end
@everywhere function kstest_transect(waitingtimes::VVVQ, r, Δt, R, Cₛ, C₀, T)
    [kstest_transect(waitingtimes[i],r[i],R,Cₛ,C₀,T) for i in eachindex(waitingtimes)]
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
        T = params["T"]u"ms"
        ks = kstest_transect(waitingtimes, r, R, Cₛ, C₀, T)
        waitingtimes = nothing
        GC.gc()
        # evaluate mean in each interval
        l = minimum(length.(ks))
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
