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

datasets_hein = unpack_dataframe(
    collect_results(datadir("HeinMod"); rinclude=[r"IC"])
)

f = jldopen(datadir("Poisson", "RC.jld2"), "r")
Rks, Cks = ustrip.(f["R"]), ustrip.(f["Cₛ"])
close(f)
datasets_ks = unpack_dataframe(
    collect_results(datadir("Poisson", "KolmogorovSmirnov"); rinclude=[r"IC"])
)

## Initialize figure layout
w, h = TwoColumns(1.5)
fig = Figure(resolution = (w,h))
pa = fig[1:2,1] = GridLayout()
pb = fig[1:2,2] = GridLayout()
pc = fig[1:2,3] = GridLayout()

## Some global constants
Rmin, Rmax = 0.1, 70.0
Cmin, Cmax = 1.5e-3, Cks[end]
cmap = :viridis
clims = (0, 3)
clevels = range(clims...; step = 0.25)

Ic_str = rich("I", subscript("c"))
logIc_str = rich("log", subscript("10"), Ic_str)

## Tall global colorbar on the right edge
cb = Colorbar(fig[1:2,4],
    colormap = cmap,#cgrad(cmap, length(clevels)-1; categorical=true),
    colorrange = clims,
    ticks = 0:4,
    ticksvisible = false,
    label = logIc_str,
    #height = 0.85*h
)

## Panel C: large Dc
D0 = 1000
U0 = 50
T0 = 100
Π0 = 6
df_hein = subset(datasets_hein,
    :Dc => Dc -> Dc .== D0,
    :U => U -> U .== U0,
    :T => T -> T .== T0,
    :Π => Π -> Π .== Π0,
)[1,:]

df_ks = subset(datasets_ks,
    :Dc => Dc -> Dc .== D0,
    :U => U -> U .== U0,
    :T => T -> T .== T0,
)[1,:]

ax_c1 = Axis(pc[1,1],
    xscale = log10,
    yscale = log10,
    xticks = [1, 3, 9, 27],
    yticks = [0.01, 0.1, 1],
    xticklabelsvisible = false,
    yticklabelsvisible = false,
    #title = rich("D", subscript("c"), " = 1000 μm²/s"),
    title = "Large Diffusivity",
    titlefont = :bold
)
xlims!(ax_c1, (Rmin, Rmax))
ylims!(ax_c1, (Cmin, Cmax))

contourf!(ax_c1,
    R, Cₛ, log10.(df_hein.ic),
    colormap = cmap,
    levels = clevels,
)
contour!(ax_c1,
    R, Cₛ, log10.(df_hein.ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 8,
    color = :white
)

ax_c2 = Axis(pc[2,1],
    xscale = log10,
    yscale = log10,
    xticks = [1, 3, 9, 27],
    yticks = [0.01, 0.1, 1],
    yticklabelsvisible = false,
    xlabel = "R (μm)"
)
xlims!(ax_c2, (Rmin, Rmax))
ylims!(ax_c2, (Cmin, Cmax))

heatmap!(ax_c2,
    Rks, Cks, log10.(df_ks.ic),
    colormap = cmap,
    colorrange = clims,
)
contour!(ax_c2,
    R, Cₛ, log10.(df_hein.ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 8,
    color = :white
)



## Panel B: large U
D0 = 500
U0 = 100
T0 = 100
Π0 = 6
df_hein = subset(datasets_hein,
    :Dc => Dc -> Dc .== D0,
    :U => U -> U .== U0,
    :T => T -> T .== T0,
    :Π => Π -> Π .== Π0,
)[1,:]

df_ks = subset(datasets_ks,
    :Dc => Dc -> Dc .== D0,
    :U => U -> U .== U0,
    :T => T -> T .== T0,
)[1,:]

ax_b1 = Axis(pb[1,1],
    xscale = log10,
    yscale = log10,
    xticks = [1, 3, 9, 27],
    yticks = [0.01, 0.1, 1],
    xticklabelsvisible = false,
    yticklabelsvisible = false,
    #title = "U = 100 μm/s",
    title = "High Speed",
    titlefont = :bold
)
xlims!(ax_b1, (Rmin, Rmax))
ylims!(ax_b1, (Cmin, Cmax))

contourf!(ax_b1,
    R, Cₛ, log10.(df_hein.ic),
    colormap = cmap,
    levels = clevels,
)
contour!(ax_b1,
    R, Cₛ, log10.(df_hein.ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 8,
    color = :white
)

ax_b2 = Axis(pb[2,1],
    xscale = log10,
    yscale = log10,
    xticks = [1, 3, 9, 27],
    yticks = [0.01, 0.1, 1],
    yticklabelsvisible = false,
    xlabel = "R (μm)"
)
xlims!(ax_b2, (Rmin, Rmax))
ylims!(ax_b2, (Cmin, Cmax))

heatmap!(ax_b2,
    Rks, Cks, log10.(df_ks.ic),
    colormap = cmap,
    colorrange = clims,
)
contour!(ax_b2,
    R, Cₛ, log10.(df_hein.ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 8,
    color = :white
)


## Panel A: small T
D0 = 500
U0 = 50
T0 = 50
Π0 = 6
df_hein = subset(datasets_hein,
    :Dc => Dc -> Dc .== D0,
    :U => U -> U .== U0,
    :T => T -> T .== T0,
    :Π => Π -> Π .== Π0,
)[1,:]

df_ks = subset(datasets_ks,
    :Dc => Dc -> Dc .== D0,
    :U => U -> U .== U0,
    :T => T -> T .== T0,
)[1,:]

ax_a1 = Axis(pa[1,1],
    xscale = log10,
    yscale = log10,
    xticks = [1, 3, 9, 27],
    yticks = [0.01, 0.1, 1],
    xticklabelsvisible = false,
    ylabel = "Cₛ (μM)",
    #title = "T = 50 ms",
    title = "Short Sensory Timescale",
    titlefont = :bold
)
xlims!(ax_a1, (Rmin, Rmax))
ylims!(ax_a1, (Cmin, Cmax))

contourf!(ax_a1,
    R, Cₛ, log10.(df_hein.ic),
    colormap = cmap,
    levels = clevels,
)
contour!(ax_a1,
    R, Cₛ, log10.(df_hein.ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 8,
    color = :white
)

ax_a2 = Axis(pa[2,1],
    xscale = log10,
    yscale = log10,
    xticks = [1, 3, 9, 27],
    yticks = [0.01, 0.1, 1],
    xlabel = "R (μm)",
    ylabel = "Cₛ (μM)"
)
xlims!(ax_a2, (Rmin, Rmax))
ylims!(ax_a2, (Cmin, Cmax))

heatmap!(ax_a2,
    Rks, Cks, log10.(df_ks.ic),
    colormap = cmap,
    colorrange = clims,
)
contour!(ax_a2,
    R, Cₛ, log10.(df_hein.ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 8,
    color = :white
)


##
fig
