using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames, Distributions, HypothesisTests

"""
    has_different_distribution(waiting_times, c, Δt)
Compare the empirical distribution of waiting times in `waiting_times` with
a reference distribution of waiting times at constant concentration `c`.
The reference distribution is a truncated exponential with lower bound `Δt`
and average `waitingtime(c)`.
"""
function has_different_distribution(waiting_times, c, Δt)
    isempty(waiting_times) && return Float64(false)
    # average waiting time at concentration c
    τ = waitingtime(c) |> u"ms" |> ustrip
    tmin = Δt |> u"ms" |> ustrip
    # expected pdf of waiting times at concentration c
    reference_pdf = Truncated(Exponential(τ), tmin, Inf)
    # compare sampled waiting times to reference pdf via Kolmogorov-Smirnov test
    KStest = ApproximateOneSampleKSTest(ustrip.(waiting_times), reference_pdf)
    # if p < 0.05 the two distributions are not the same
    return Float64(pvalue(KStest) < 0.05)
end

function sensing(R, L, params)
    @unpack C₀, T, U, Δt, N = params
    iter = Iterators.product(eachindex(R), eachindex(L)) |> collect
    diffsensing = [Float64[] for _ in iter]
    grids = [range(r,r,length=1) for r in R]
    for itr in iter
        i,j = itr # R, L
        in_params = @strdict C₀ T U Δt N i j
        in_params = Dict(keys(in_params) .=> ustrip.(values(in_params)))
        dataset = jldopen(datadir("PoissonSampling", savename("waitingtimes", in_params, "jld2")))
        r = dataset["r"]
        waiting_times = dataset["waiting_times"]
        diffsensing[i,j] = zeros(length(r))
        grids[i] = r
        for k in eachindex(r)
            # concentration at the midpoint of current interval
            c = C(r[k],L[j],C₀)
            # fraction of samples where waiting times dont match the reference
            ϕ = mean(has_different_distribution.(waiting_times[k], c, Δt))
            diffsensing[i,j][k] = ϕ
        end
    end
    grids, diffsensing
end

function getksradius(grid, diffsensing, R)
    idx = findfirst(diffsensing .> 0.5)
    isnothing(idx) ? R : grid[idx]
end

function phycosphere_ks()
    f = jldopen(datadir("PoissonSampling", "paramspaceRL.jld2"))
    R, L = f["R"], f["L"]
    # parameters for the sensing process
    C₀ = 1u"nM" # background concentration
    T = 100u"ms" # sensory integration timescale
    U = 46.5u"μm/s" # speed
    Δt = 1e-4u"ms" # temporal resolution of the simulation
    N = 100
    unit_params = @strdict C₀ T U Δt N
    params = Dict(keys(unit_params) .=> ustrip.(values(unit_params)))

    # spatially-resolved sensing statistics
    sensing_data = produce_or_load(
        datadir("PoissonSampling", "KolmogorovSmirnov"), params;
        prefix="sensing", suffix="jld2", tag=false
    ) do params
        grids, diffsensing = sensing(R, L, unit_params)
        @strdict grids diffsensing
    end
    
    # phycosphere estimation
    produce_or_load(
        datadir("PoissonSampling", "KolmogorovSmirnov"), params;
        prefix="phycosphere", suffix="jld2", tag=false
    ) do params
        # [1] takes the data, [2] is path
        @unpack grids, diffsensing = sensing_data[1]
        S = [
            getksradius(grids[i], diffsensing[i,j], R[i])
            for i in eachindex(R), j in eachindex(L)
        ]
        @strdict S
    end
end

phycosphere_ks()