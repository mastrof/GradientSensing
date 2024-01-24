##
using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using GLMakie
using PublicationFiguresMakie
arial_italic = "/usr/share/fonts/TTF/ariali.ttf"
set_theme!(Publication,
    Axis = (
        xminorticksvisible=false, yminorticksvisible=false,
        xtickcolor = :white, ytickcolor = :white
    ),
    fonts = (
        regular = "Arial",
        bold = "Arial Bold",
        italic = arial_italic,
    ) 
)

function wsmooth(x)
    w = [1 1 1; 1 4 1; 1 1 1]
    replacenans(y) = [isnan(z) ? zero(z) : z for z in y]
    mapwindow(y -> mean(replacenans(y) .* w), x, (3,3))
end

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
Rmin, Rmax = 0.36, 45.0
Cmin, Cmax = 1.5e-3, Cks[end]
cmap = :viridis
clims = (0, 2.5)
clevels = range(clims...; step = 0.25)

Ic_str = rich("I", subscript("c"); font=:italic)
C_str = rich(rich("C", subscript("s"); font=:italic), " (μM)")
R_str = rich(rich("R"; font=:italic), " (μm)")

df_ref = subset(datasets_hein,
    :Dc => Dc -> Dc .== 500,
    :U => U -> U .== 50,
    :T => T -> T .== 100,
    :Π => Π -> Π .== 6
)[1,:]

## Tall global colorbar on the right edge
cb = Colorbar(fig[1:2,4],
    colormap = cmap,#cgrad(cmap, length(clevels)-1; categorical=true),
    colorrange = clims,
    ticks = 0:4,
    ticksvisible = false,
    tickformat = values -> [rich("10", superscript("$(Int(z))")) for z in values],
    label = Ic_str,
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
    xlabel = R_str
)
xlims!(ax_c2, (Rmin, Rmax))
ylims!(ax_c2, (Cmin, Cmax))

contourf!(ax_c2,
    Rks, Cks, log10.(wsmooth(df_ks.ic)),
    colormap = cmap,
    colorrange = clims,
    levels = clevels
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
    xlabel = R_str
)
xlims!(ax_b2, (Rmin, Rmax))
ylims!(ax_b2, (Cmin, Cmax))

contourf!(ax_b2,
    Rks, Cks, log10.(wsmooth(df_ks.ic)),
    colormap = cmap,
    colorrange = clims,
    levels = clevels
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
    ylabel = C_str,
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
    xlabel = R_str,
    ylabel = C_str
)
xlims!(ax_a2, (Rmin, Rmax))
ylims!(ax_a2, (Cmin, Cmax))

contourf!(ax_a2,
    Rks, Cks, log10.(wsmooth(df_ks.ic)),
    colormap = cmap,
    colorrange = clims,
    levels = clevels
)
contour!(ax_a2,
    R, Cₛ, log10.(df_hein.ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 8,
    color = :white
)

## Reference line in all panels
for ax in (ax_a1, ax_a2, ax_b1, ax_b2, ax_c1, ax_c2)
    contour!(ax,
    R, Cₛ, log10.(df_ref.ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 5,
    color = RGBAf(1.0, 0.75, 0.8, 0.4)
)

end

## Slightly increase column gap
colgap!(fig.layout, Relative(0.03))

##
fig
