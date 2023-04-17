using DrWatson
@quickactivate :GradientSensing
using JLD2

function spatial_counting(x₀, x₁, L, C₀, U, Δt)
    event_times = zero([Δt])
    # use expected counts at midpoint to sizehint event_times
    xp = (x₀+x₁)/2
    n = nevents(C(xp,L,C₀), Δt)
    sizehint!(event_times, 2*ceil(Int,n))
    # start from x₀ and move in steps U*Δt
    i = 0
    x = x₀
    while x < x₁
        i += 1
        x += U*Δt
        ω = eventrate(C(x,L,C₀))
        if rand() < upreferred(ω*Δt)
            # if an event occurs, record its timing
            push!(event_times, i*Δt)
        end
    end
    diff(event_times)
end

function sample_waitingtimes(r, L, C₀, U, Δt, N)
    # initialize nested arrays to hold all samples
    waiting_times = [[zero([Δt]) for _ in 1:N] for _ in eachindex(r)[2:end]]
    # sample poisson events in each interval of the landscape
        # and repeat N times for each interval
        for k in eachindex(r)[2:end]
            for n in 1:N
            # !! remember r[k] < r[k-1] !!
            x₀, x₁ = r[k], r[k-1]
            waiting_times[k-1][n] = spatial_counting(x₀, x₁, L, C₀, U, Δt)
        end
    end
    waiting_times
end

function poisson_sampling()
    f = jldopen(datadir("PoissonSampling", "paramspaceRL.jld2"))
    R, L = f["R"], f["L"]
    # parameters for the sensing process
    C₀ = 1u"nM" # background concentration
    T = 100u"ms" # sensory integration timescale
    U = 46.5u"μm/s" # speed
    Δx = U*T |> u"μm" # spatial resolution of the sensor
    Δt = 1e-4u"ms" # temporal resolution of the simulation
    N = 100 # number of repetitions for each sampling

    # for each R define a grid from 20R to R in steps Δx
    # !! notice this starts away from the source and moves towards it !!
    rs = [reverse(range(Rᵢ, max(100u"μm",20Rᵢ), step=Δx)) for Rᵢ in R]

    iter = Iterators.product(eachindex(R), eachindex(L)) |> collect
    Threads.@threads for itr in iter
        i, j = itr # R, L
        # initialize nested arrays to hold all samples
        waiting_times = sample_waitingtimes(rs[i], L[j], C₀, U, Δt, N)

        params = @strdict C₀ T U Δt N i j
        # remove units from parameters to use for savename
        params = Dict(keys(params) .=> ustrip.(values(params)))
        produce_or_load(
            datadir("PoissonSampling"), params;
            prefix="waitingtimes", suffix="jld2",
            tag=false
        ) do params
            local r = midpoints(rs[i])
            @strdict waiting_times r
        end
    end
end

poisson_sampling()