##
using Distributed
@everywhere using DrWatson
@everywhere @quickactivate :GradientSensing
@everywhere begin
    using JLD2, DataFrames, Unitful
end

## aliases for brevity
@everywhere begin
    VVQ = AbstractVector{<:AbstractVector{<:Quantity}}
    VVVQ = AbstractVector{<:AbstractVector{<:AbstractVector{<:Quantity}}}
end

##
@everywhere function produce_data()
    datasets = get_filelist_waitingtimes()
    pmap(gradient_estimation, datasets)
    GC.gc()
end

@everywhere function get_filelist_waitingtimes()
    filenames = readdir(
        datadir("Poisson"); # on local
        #joinpath(ENV["SCRATCH"], "Poisson"); # on cluster
        join = true
    )
    filter!(s -> contains(s, "waitingtimes"), filenames)
    filenames
end

@everywhere function gradient_estimation(fname)
    params = parse_savename(fname)[2]
    produce_or_load(
        datadir("Poisson", "LinearRegression"), params; # on local
        #joinpath(ENV["SCRATCH"], "Poisson", "LinearRegression"), params; # on cluster
        prefix="sensing", suffix="jld2", tag=false, loadfile=false
    ) do params
        f = jldopen(fname, "r")
        r = f["r"]
        waitingtimes = f["waitingtimes"]
        close(f)
        Dc = params["Dc"]u"μm^2/s"
        T = params["T"]u"ms"
        N = params["N"]u"1"
        grid, concentration, gradient = sensing(waitingtimes, r, T, Dc, N)
        @strdict grid concentration gradient
    end
    return nothing
end

@everywhere function sensing(waitingtimes::VVVQ, r, T, Dc, N)
    l = minimum(length.(r))
    grid = zeros(typeof(1.0u"μm"), l)
    concentration = zeros(typeof(1.0u"μM"), l)
    gradient = zeros(typeof(1.0u"μM/s"), l)
    for i in eachindex(waitingtimes)
        _r, _c, _g = sensing(waitingtimes[i], r[i], T, Dc)
        grid .+= _r[end-l+1:end]
        concentration .+= _c[end-l+1:end]
        gradient .+= _g[end-l+1:end]
    end
    grid ./= N
    concentration ./= N
    gradient ./= N
    grid, concentration, gradient
end

@everywhere function sensing(waitingtimes::VVQ, r, T, Dc)
    concentration = zeros(typeof(1.0u"μM"), length(r))
    gradient = zeros(typeof(1.0u"μM/s"), length(r))
    for k in eachindex(r)
        n = length(waitingtimes[k])
        concentration[k] = n / (4π*Dc*a*T*Unitful.Na) |> u"μM"
        abstimes = cumsum(waitingtimes[k]) # event times in [0,T]
        t_shifted = abstimes .- (T/2) # event times in [-T/2,+T/2]
        t = sum(t_shifted)
        gradient[k] = -12*t / (4π*Dc*a*T^3*Unitful.Na) |> u"μM/s"
    end
    r, concentration, gradient
end

##
produce_data()
