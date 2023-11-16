##
using DrWatson
@quickactivate :GradientSensing
using DelimitedFiles, JLD2, DataFrames
using Random
using GLMakie
using PublicationFiguresMakie
set_theme!(Publication,
    Axis = (
        xminorticksvisible = false,
        yminorticksvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        xticklabelsvisible = false,
    )
)

##
function generate_grid(R, dx, rng)
    d = (rand(rng) - 1/2) * dx
    x_end = R + d
    transect_length = max(100u"μm", 20R)
    x_start = x_end + transect_length
    return reverse(range(x_end, x_start, step=dx))
end

# r[1] > r[end] !!
function spatial_counting(r, R, Cs, C0, Dc, U, dt, rng)
    map(
        k -> spatial_counting(r[k], r[k-1], R, Cs, C0, Dc, U, dt, rng),
        eachindex(r)[2:end]
    )
end

function spatial_counting(x0, x1, R, Cs, C0, Dc, U, dt, rng)
    event_times = zero([dt])
    # use expected counts at midpoint to sizehint event_times
    xp = (x0 + x1) / 2
    n = nevents(C(xp, R, Cs, C0), dt, Dc)
    sizehint!(event_times, 2*ceil(Int, n))
    # start from x0 and move in steps U*dt
    i = 0
    x = x0
    while x < x1
        i += 1
        x += U*dt
        y = abs(x) < R ? R : abs(x) # flat profile inside source
        w = eventrate(C(y, R, Cs, C0), Dc)
        @assert ustrip(w) > 0
        if rand(rng) < upreferred(w*dt)
            push!(event_times, i*dt)
        end
    end
    return diff(event_times)
end

## small source
rng = Xoshiro(3)
R = 1.0u"μm"
U = 50u"μm/s"
T = 100u"ms"
C0 = 0.001u"μM"
Cs = 0.0025u"μM"
Dc = 500.0u"μm^2/s"
dt = 1e-4u"ms"

offset = U*T/1.3 # some arbitrary offset to shift the grid across the source
grid = generate_grid(R, U*T, rng) .- offset .|> u"μm"
waitingtimes = spatial_counting(grid, R, Cs, C0, Dc, U, dt, rng) .|> ustrip
nwindows = 7
z = [reverse(waitingtimes[end-i]) for i in nwindows-1:-1:0]
tevents = cumsum(vcat(z...))

fig = Figure(resolution = (721.444, 214.819))
ax = Axis(fig[1,1])
subsampling = 20
vlines!(ax, tevents[1:subsampling:end])
x1, x2 = round.(Int, extrema(tevents))
dx = (x2-x1) / nwindows
for j in 2:nwindows
    vlines!(ax, (j-1)*dx, color=:orange, linewidth=7, linestyle=:dash)
end
tightlimits!(ax)
fig
save(plotsdir("tevents_small.svg"), fig)

## large source
rng = Xoshiro(3)
R = 20.0u"μm"
U = 50u"μm/s"
T = 100u"ms"
C0 = 0.001u"μM"
Cs = 0.01u"μM"
Dc = 500.0u"μm^2/s"
dt = 1e-4u"ms"

offset = -5u"μm" # some arbitrary offset to shift the grid away from the source
grid = generate_grid(R, U*T, rng) .- offset .|> u"μm"
waitingtimes = spatial_counting(grid, R, Cs, C0, Dc, U, dt, rng) .|> ustrip
nwindows = 7
z = [reverse(waitingtimes[end-i]) for i in nwindows-1:-1:0]
tevents = cumsum(vcat(z...))

fig = Figure(resolution = (721.444, 214.819))
ax = Axis(fig[1,1])
subsampling = 50
vlines!(ax, tevents[1:subsampling:end])
x1, x2 = round.(Int, extrema(tevents))
dx = (x2-x1) / nwindows
for j in 2:nwindows
    vlines!(ax, (j-1)*dx, color=:orange, linewidth=7, linestyle=:dash)
end
tightlimits!(ax)
fig
save(plotsdir("tevents_large.svg"), fig)
