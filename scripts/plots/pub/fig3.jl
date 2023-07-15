##
using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using CairoMakie
using LaTeXStrings
using PublicationFiguresMakie
set_theme!(Publication,
    Axis=(
        xminorticksvisible=false, yminorticksvisible=false,
    )
)


## Load data
f = jldopen(datadir("HeinMod", "RC.jld2"), "r")
R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
close(f)

alldatasets = unpack_dataframe(
    collect_results(datadir("HeinMod"), rinclude=[r"IC"])
)

D̄ = 500
Ū = 50
T̄ = 100
Π̄ = 6
df = subset(alldatasets,
    :Dc => Dc -> Dc .== D̄,
    :U => U -> U .== Ū,
    :T => T -> T .== T̄,
    :Π => Π -> Π .== Π̄,
)[1,:]

## Collect data for cuts along R and C
js = map(r -> findfirst(R .≥ r), [1.0, 2.5, 5.0, 10.0, 25.0])
Rx = R[js]
I_R = df.ic[js,:]

js = map(c -> findfirst(Cₛ .≥ c), [0.015, 0.03, 0.05, 0.1])
Cx = Cₛ[js]
I_C = df.ic[:,js]'

## Some global constants
Rmin, Rmax = 0.5, 70.0
Cmin, Cmax = 2e-3, 1.0

## Rich-text strings
Ic_str = rich("I", subscript("c"))
logIc_str = rich("log", subscript("10"), Ic_str)

## Initialize figure layout
fig = Figure(resolution = TwoColumns(1.25))
pa = fig[1:2,1] = GridLayout()
pb = fig[1,2] = GridLayout()
pc = fig[2,2] = GridLayout()


## Panel B -- iso-R cuts
cmap = palette(:thermal, length(Rx))
labels = string.(round.(Rx, sigdigits=2))

ax_b = Axis(pb[1,1],
    xscale = log10,
    yscale = log10,
    xlabel = "Cₛ (μM)",
    ylabel = Ic_str,
    xticks = [0.01, 0.1, 1],
    yticks = [1, 10, 100]
)
series!(ax_b,
    Cₛ, I_R,
    color = cmap,
    linewidth = 9,
    labels = labels,
)
axislegend("R (μm)"; position=:lt)
xlims!(ax_b, (Cmin, Cₛ[end]))
ylims!(ax_b, (1, 100))


## Panel C -- iso-C cuts
cmap = palette(:batlow, length(Cx))
labels = string.(round.(Cx, sigdigits=2))

ax_c = Axis(pc[1,1],
    xscale = log10,
    yscale = log10,
    xlabel = "R (μm)",
    ylabel = Ic_str,
    xticks = [1, 3, 9, 27],
    yticks = [1, 10, 100]
)
series!(ax_c,
    R, I_C,
    color = cmap,
    linewidth = 9,
    labels = labels
)
axislegend("Cₛ (μM)"; position=:rt)
xlims!(ax_c, (Rmin, Rmax))
ylims!(ax_c, (1, 100))


## Panel A -- IC landscape
ax_a = Axis(pa[1,1],
    xlabel = "R (μm)",
    ylabel = "Cₛ (μM)",
    xscale = log10,
    yscale = log10,
    xticks = [1, 3, 9, 27],
    yticks = [0.01, 0.1, 1]
)

clims = (0, 3)
clevels = range(clims...; step = 0.25)
cmap = :viridis

cb = Colorbar(pa[1,2],
    colormap = cgrad(cmap, length(clevels)-1; categorical=true),
    colorrange = clims,
    ticks = 0:4,
    ticksvisible = false,
    label = logIc_str
)

contourf!(ax_a,
    R, Cₛ, log10.(df.ic),
    colormap = cmap,
    levels = clevels,
    extendlow = :white
)
contour!(ax_a,
    R, Cₛ, log10.(df.ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 8,
    color = :white
)

scatter!(ax_a,
    Rx, repeat([1*Cmin], length(Rx)),
    marker = :utriangle,
    markersize = 35,
    strokewidth = 1,
    strokecolor = :white,
    color = palette(:thermal, length(Rx))
)

scatter!(ax_a,
    repeat([1*Rmin], length(Cx)), Cx,
    marker = :rtriangle,
    markersize = 35,
    strokewidth = 1,
    strokecolor = :white,
    color = palette(:batlow, length(Cx))
)

xlims!(ax_a, (Rmin, Rmax))
ylims!(ax_a, (Cmin, Cmax))

##
fig
