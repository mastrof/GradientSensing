##
using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using CairoMakie
using PublicationFiguresMakie
set_theme!(Publication,
    Axis=(
        xminorticksvisible=false, yminorticksvisible=false,
        xtickcolor = :white, ytickcolor=:white,
    )
)


## Load data
f = jldopen(datadir("HeinMod", "RC.jld2"), "r")
R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
close(f)

Ū = 50
T̄ = 100
D̄ = 500
Π̄ = 6

alldatasets = unpack_dataframe(
    collect_results(datadir("HeinMod"), rinclude=[r"IC"])
)
df = subset(alldatasets,
    :U => U -> U .== Ū,
    :T => T -> T .== T̄,
    :Dc => Dc -> Dc .== D̄,
    :Π => Π -> Π .== Π̄
)[1,:]

## Phytoplankton scaling
PER = [0.02, 0.1, 0.5]
r = R[0.5 .< R .< 60]
L = leakage_rate.(r .* 1u"μm", PER')
Cp = C.(r .* 1u"μm", L, 1u"nM", D̄*1u"μm^2/s")
Cp = ustrip.(Cp)

## Initialize figure layout
w, h = SinglePlot()
fig = Figure(resolution = (w,h))

## Some global constants
Rmin, Rmax = r[1], r[end]
Cmin, Cmax = 1.5e-3, 1.0
cmap = :viridis
clims = (0, 3)
clevels = range(clims...; step = 0.25)

Ic_str = rich("I", subscript("c"))
logIc_str = rich("log", subscript("10"), Ic_str)

## Tall global colorbar on the right edge
cb = Colorbar(fig[1,2],
    colormap = cgrad(cmap, length(clevels)-1; categorical=true),
    colorrange = clims,
    ticks = 0:4,
    ticksvisible = false,
    label = logIc_str,
)

## Panel
ax = Axis(fig[1,1],
    xscale = log10,
    yscale = log10,
    xticks = [1, 3, 9, 27],
    yticks = [0.01, 0.1, 1],
    xlabel = "R (μm)",
    ylabel = "Cₛ (μM)"
)
ylims!(Cmin, Cmax)
xlims!(Rmin, Rmax)

ic = df.ic
contourf!(ax, R, Cₛ, log10.(ic),
    colormap = cmap,
    levels = clevels,
)

contour!(ax, R, Cₛ, log10.(df.ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 8,
    color = :white
)

band!(ax, r, Cp[:,1], Cp[:,3], color=RGBAf(1,0.5,0.7,0.25))
lines!(ax, r, Cp[:,2], linewidth=8, color=RGBAf(1,0.5,0.7,0.95))


##
fig
