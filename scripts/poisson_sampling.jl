using DrWatson
@quickactivate :GradientSensing
using JLD2, Base.Threads, ChunkSplitters

"""
    spatial_counting(x₀, x₁, R, Cₛ, C₀, U, Δt)
Simulates a concentration sensor moving with speed `U` in an interval [x₀,x₁]
while counting Poissonian absorption events from the surrounding concentration field.
The concentration field is produced by a spherical source with radius `R`
and excess surface concentration `Cₛ` whose center is located at the origin (x=0).

Returns a vector of waiting times between successive counting events.

**Arguments**
- `x₀`: initial position
- `x₁`: final position (`x₁>x₀`)
- `R`: source radius
- `Cₛ`: excess concentration at the surface of the source
- `C₀`: backgrund concentration (at infinite distance from the source)
- `U`: velocity with which the sensor moves through the interval
- `Δt`: minimal temporal resolution at which the sensor can distinguish two signals
"""
function spatial_counting(x₀, x₁, R, Cₛ, C₀, U, Δt)
    event_times = zero([Δt])
    # use expected counts at midpoint to sizehint event_times
    xp = (x₀+x₁)/2
    n = nevents(C(xp,R,Cₛ,C₀), Δt)
    sizehint!(event_times, 2*ceil(Int,n))
    # start from x₀ and move in steps U*Δt
    i = 0
    x = x₀
    while x < x₁
        i += 1
        x += U*Δt
        ω = eventrate(C(x,R,Cₛ,C₀))
        if rand() < upreferred(ω*Δt)
            # if an event occurs, record its timing
            push!(event_times, i*Δt)
        end
    end
    diff(event_times)
end

"""
    spatial_counting(r, R, Cₛ, C₀, U, Δt)
Simulates a concentration sensor moving with speed `U` along a grid `r`.
while counting Poissonian absorption events from the surrounding concentration field.
The concentration field is produced by a spherical source with radius `R`
and excess surface concentration `Cₛ` whose center is located at the origin (x=0).

Returns a vector of waiting times between successive counting events for each
interval that makes up the grid `r`.
Measurements are independent between intervals.

**Arguments**
- `r: spatial grid with mesh size equal to the spatial resolution of the sensor; must be sorted in reverse order (`r[1]>r[end]`) 
- `R`: source radius
- `Cₛ`: excess concentration at the surface of the source
- `C₀`: backgrund concentration (at infinite distance from the source)
- `U`: velocity with which the sensor moves through the interval
- `Δt`: minimal temporal resolution at which the sensor can distinguish two signals
"""
function spatial_counting(r, R, Cₛ, C₀, U, Δt)
    map(k -> spatial_counting(r[k], r[k-1], R, Cₛ, C₀, U, Δt), eachindex(r)[2:end])
end

"""
    sample_waitingtimes(params)
Performs spatial sampling of a concentration profile by measuring waiting times
between Poissonian absorption events along a transect towards a spherical source.

The `params` dictionary must have the following keys, which define the
properties of the simulation:
- `R`: source radius
- `Cₛ`: excess concentration at the surface of the source
- `C₀`: background concentration (at infinite distance from the source)
- `U`: speed at which the sensor moves along the transect
- `T`: sensory integration timescale of the sensor
- `Δt`: minimal temporal resolution at which the sensor can distinguish two events
- `N`: number of repetitions for each sampling transect

For any given set of parameters, a grid is generated which starts away from
the source and moves towards it in steps of size `Δx=U*T`.
Within any spatial interval `Δx`, the sensor resolves a number of absorption
events which is determined by the local concentration field arising from the source.
Measurements between intervals are independent (i.e. the state of the sensor is
reset in every spatial interval).
Each repetition is stored separately.

The data is saved to a jld2 file containing the objects `waitingtimes` and `r`.

`r`: the midpoints of the spatial grid which has been sampled by the sensor

`waitingtimes`: a `Vector{Vector{Vector{Quantity{<:Real,Unitful.𝐓}}}}`
containing the waiting times sampled during experiments.
`waitingtimes[i][j]` returns the list of waiting times associated to the `i`-th
repetition of the experiment at position `r[j]`.
"""
function sample_waitingtimes(params)
    # remove units for savename
    params_ustrip = Dict(keys(params) .=> ustrip.(values(params)))
    produce_or_load(
        datadir("PoissonSampling"), params_ustrip;
        prefix="waitingtimes", suffix="jld2",
        tag=false
    ) do params_ustrip
        @unpack R, Cₛ, C₀, U, T, Δt, N = params

        Δx = U*T |> u"μm" # spatial resolution of the sensor
        # define a grid from max(100u"μm",20R) to R in steps Δx
        # !! notice this starts away from the source and moves towards it !!
        r = reverse(range(R, max(100u"μm",20R), step=Δx))
        waitingtimes = map(_ -> spatial_counting(r, R, Cₛ, C₀, U, Δt), 1:N)
        r = midpoints(r)
        @strdict waitingtimes r
    end
end

"""
Performs a molecule counting experiment for each `(R,Cₛ)` pair in
the parameter space. (See `sample_waitingtimes`.)
"""
function poisson_sampling()
    f = jldopen(datadir("PoissonSampling", "paramspaceRC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]
    # parameters for the sensing process
    T = 100u"ms" # sensory integration timescale
    U = 46.5u"μm/s" # speed
    Δt = 1e-6u"ms" # temporal resolution of the simulation
    N = 100 # number of repetitions for each sampling

    iter = Iterators.product(R, Cₛ) |> collect
    @sync for (k_range, _) in chunks(iter, nthreads(), :scatter)
        @spawn for k in k_range
            local R, Cₛ = iter[k]
            params = @strdict R Cₛ C₀ T U Δt N
            sample_waitingtimes(params)
        end
    end
end

poisson_sampling()