using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames

function sensing(R, Cₛ, params)
    @unpack C₀, T, U, Δt, N = params
    iter = Iterators.product(eachindex(R), eachindex(Cₛ)) |> collect
    concentration = [empty([Float64(C₀)]) for _ in iter]
    gradient = [empty([Float64(C₀/T)]) for _ in iter]
    grids = [range(r,r,length=1) for r in R]
    for itr in iter
        i,j = itr # R, Cₛ
        in_params = @strdict C₀ T U Δt N i j
        in_params = Dict(keys(in_params) .=> ustrip.(values(in_params)))
        dataset = jldopen(datadir("Poisson", savename("waitingtimes", in_params, "jld2")))
        r = dataset["r"]
        waiting_times = dataset["waiting_times"]
        concentration[i,j] = zeros(typeof(C₀), length(r))
        gradient[i,j] = zeros(typeof(C₀/T), length(r))
        grids[i] = r
        for k in eachindex(r)
            # estimated concentration (Mora-Wingreen 2010 eq. S29)
            n = mean(length.(waiting_times[k]))
            concentration[i,j][k] = n/(4π*Dc*a*T*Unitful.Na) |> u"μM"
            # estimated (temporal) gradient (Mora-Wingreen 2010 eq. S30)
            # event times should be shifted in the [-T/2,T/2] interval
            abstimes = cumsum.(waiting_times[k]) # event times in [0,T]
            t_shifted = map(ts -> (T/2).-ts, abstimes) # event times in [-T/2,T/2]
            t = mean(sum.(t_shifted))
            gradient[i,j][k] = 12*t/(4π*Dc*a*T^3*Unitful.Na) |> u"μM/s"
        end
    end
    grids, concentration, gradient
end

function gradient_estimation()
    f = jldopen(datadir("Poisson", "RL.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]
    # parameters for the sensing process
    C₀ = 1u"nM" # background concentration
    T = 100u"ms" # sensory integration timescale
    U = 46.5u"μm/s" # speed
    Δt = 1e-4u"ms" # temporal resolution of the simulation
    N = 100
    unit_params = @strdict C₀ T U Δt N
    params = Dict(keys(unit_params) .=> ustrip.(values(unit_params)))

    sensing_data = produce_or_load(
        datadir("Poisson", "LinearRegression"), params;
        prefix="sensing", suffix="jld2", tag=false
    ) do params
        grids, concentration, gradient = sensing(R, L, unit_params)
        @strdict grids concentration gradient
    end
end

gradient_estimation()
