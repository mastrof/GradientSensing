export setup_abm, KSMicrobe, ks_sensing!

function setup_abm(n::Int, R::Real, Cₛ::Real, C₀::Real, T::Real, Δt::Real, U::Real, L::Real, D::Int=3)
    space = ContinuousSpace(ntuple(_ -> L, D); periodic=true)
    properties = Dict(
        :C₀ => C₀, :Cₛ => Cₛ, :R => R,
        :concentration_field => (_,_) -> Cₛ,
    )
    model = StandardABM(KSMicrobe{D}, space, Δt; properties)
    sensory_timescale = T
    motility = RunTumble(speed=[U])
    foreach(_ -> add_agent!(model; sensory_timescale, motility), 1:n)
    model → ks_sensing!
    model
end

mutable struct KSMicrobe{D} <: AbstractMicrobe{D}
    id::Int
    pos::NTuple{D,Float64}
    motility::AbstractMotility
    vel::NTuple{D,Float64}
    speed::Float64
    turn_rate::Float64
    rotational_diffusivity::Float64
    radius::Float64
    state::Int
    sensory_timescale::Float64
    sensor::Vector{Float64}
    sensor_history::Vector{Float64}

    KSMicrobe{D}(
        id::Int = rand(1:typemax(Int32)),
        pos::NTuple{D,<:Real} = ntuple(zero, D);
        motility::AbstractMotility = RunTumble(speed=[46.5]),
        vel::NTuple{D,<:Real} = rand_vel(D),
        speed::Real = rand_speed(motility),
        turn_rate::Real = 0.5,
        rotational_diffusivity::Real = 0.0,
        radius::Real = 0.5,
        state::Int = 0,
        sensory_timescale::Real = 100,
        sensor::Vector{<:Real} = Float64[],
        sensor_history::Vector{<:Real} = Float64[]
    ) where {D} = new{D}(
        id, Float64.(pos), motility, Float64.(vel), Float64(speed),
        Float64(turn_rate), Float64(rotational_diffusivity),
        Float64(radius), state, Float64(sensory_timescale),
        Float64.(sensor), Float64.(sensor_history)
    )
end

#== Used in microbe_step! ==#
function MicrobeAgents._affect!(microbe::KSMicrobe, model)
    absorption_event!(microbe, model)
    microbe.state += 1
end

"""
Measures local concentration, estimates rate of absorption events
and if an absorption event occurs stores its time into
`microbe.sensor`.
"""
function absorption_event!(microbe::KSMicrobe, model)
    Δt = model.timestep
    c = model.concentration_field(microbe.pos, model)
    e = eventrate(c*1u"μM") |> u"s^-1" |> ustrip
    if rand(abmrng(model)) < e*Δt
        push!(microbe.sensor, model.t*Δt)
    end
end

#== Used in model.update! ==#
"""
Performs a one-tailed KS test every time interval T,
comparing the current timeseries of absorption events (`microbe.sensor`)
with the previous one (`microbe.sensor_history`).
If the new timeseries has cdf larger than the old one, the test is passed
and the microbe is removed from the model.
"""
function ks_sensing!(model)
    Δt = model.timestep
    for microbe in allagents(model)
        T = microbe.sensory_timescale
        microbe.state*Δt ≠ T && continue
        KStest = ApproximateTwoSampleKSTest(microbe.sensor, microbe.sensor_history)
        if pvalue(KStest, tail=:right) < 0.05
            remove_agent!(microbe, model)
            continue
        end
        reset_sensing!(microbe, model)
    end
end

"""
Resets the sensing state of `microbe`.
"""
function reset_sensing!(microbe::KSMicrobe, model)
    microbe.sensor_history = copy(microbe.sensor)
    microbe.sensor = Float64[]
    sizehint!(microbe.sensor, microbe.sensor_history)
    microbe.state = 0
end
