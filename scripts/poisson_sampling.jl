using Distributed
@everywhere using DrWatson
@everywhere @quickactivate :GradientSensing
@everywhere using JLD2


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

@everywhere function sample_waitingtimes(params)
    # remove units for savename
    params_ustrip = Dict(keys(params) .=> ustrip.(values(params)))
    produce_or_load(
        datadir("Poisson"), # this is valid on local
        # joinpath(ENV["SCRATCH"], "Poisson"), # this is valid on cluster
        params_ustrip;
        prefix="waitingtimes", suffix="jld2",
        tag=false, loadfile=false
    ) do params_ustrip
        produce_data(params)
    end
    GC.gc()
    return nothing
end

@everywhere function produce_data(params)
    @unpack R, Cₛ, C₀, U, T, Δt, N = params
    Δx = U*T |> u"μm"
    staggered_grids = [generate_grid(R, Δx) for _ in 1:N]
    waitingtimes = map(grid -> spatial_counting(grid, R, Cₛ, C₀, U, Δt), staggered_grids)
    r = map(midpoints, staggered_grids)
    @strdict waitingtimes r
end

"""
    generate_grid(R, Δx)
Generate a regular grid of spacing `Δx` starting at
`max(100u"μm", 20R) ± δ` and ending at `R ± δ`,
where `δ` is a random value in `(-Δx/2,+Δx/2)`.
"""
@everywhere function generate_grid(R, Δx)
    x_end = R + rand(Uniform(-1/2,+1/2))*Δx
    transect_length = max(100u"μm", 20R)
    x_start = x_end + transect_length
    return range(x_start, x_end, -Δx)
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

@everywhere function spatial_counting(r, R, Cₛ, C₀, U, Δt)
    map(k -> spatial_counting(r[k], r[k-1], R, Cₛ, C₀, U, Δt), eachindex(r)[2:end])
end

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

@everywhere function spatial_counting(x₀, x₁, R, Cₛ, C₀, U, Δt)
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
        y = x < R ? R : x # equivalent to getting stuck at the surface
        ω = eventrate(C(y,R,Cₛ,C₀))
        if rand() < upreferred(ω*Δt)
            # if an event occurs, record its timing
            push!(event_times, i*Δt)
        end
    end
    diff(event_times)
end

f = jldopen(datadir("Poisson", "RC.jld2"))
R, Cₛ = f["R"], f["Cₛ"]
# parameters for the sensing process
T = [100]u"ms" # sensory integration timescale
U = [46.5]u"μm/s" # speed
Δt = [1e-4]u"ms" # temporal resolution of the simulation
N = [100] # number of repetitions for each sampling

allparams = @strdict R Cₛ T U Δt N
dicts = dict_list(allparams)
foreach(dict -> dict["C₀"] = C₀, dicts)
pmap(sample_waitingtimes, dicts)
