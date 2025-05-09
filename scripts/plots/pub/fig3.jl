##
using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using GLMakie
using ColorSchemes
palette(colors::Symbol, n::Int) = get(colorschemes[colors], range(0,1;length=n))

## boundaries
leftboundary(R,C0,U,T,Dc,a) = 34.4*U/(a*Dc*R*Unitful.Na)
rightboundary(R,C0,U,T,Dc,a) = 34.4*R^2/(a*Dc*U^2*T^3*Unitful.Na)

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
Ic_str = rich("I", subscript("c"); font=:italic)
C_str = rich(rich("C", subscript("s"); font=:italic), " (μM)")
R_str = rich(rich("R"; font=:italic), " (μm)")


## Initialize figure layout
fig = Figure(resolution = (1600,750))
pa = fig[1:2,1] = GridLayout()
pb = fig[1,2] = GridLayout()
pc = fig[2,2] = GridLayout()


## Panel B -- iso-R cuts
cmap = palette(:thermal, length(Rx))
labels = string.(round.(Rx, sigdigits=2))

ax_b = Axis(pb[1,1],
    xscale = log10,
    yscale = log10,
    xlabel = C_str,
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
axislegend(ax_b; position=(0.02, 0))
text!(ax_b, 0.04, 0.8; text=R_str, space=:relative)
xlims!(ax_b, (Cmin, Cₛ[end]))
ylims!(ax_b, (1, 100))


## Panel C -- iso-C cuts
cmap = palette(:batlow, length(Cx))
labels = string.(round.(Cx, sigdigits=2))

ax_c = Axis(pc[1,1],
    xscale = log10,
    yscale = log10,
    xlabel = R_str,
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
axislegend(ax_c; position=(1, 0.2))
text!(ax_c, 0.8, 0.72; text=C_str, space=:relative)
xlims!(ax_c, (Rmin, Rmax))
ylims!(ax_c, (1, 100))


## Panel A -- IC landscape
ax_a = Axis(pa[1,1],
    xlabel = R_str,
    ylabel = C_str,
    xscale = log10,
    yscale = log10,
    xticks = [1, 3, 9, 27],
    yticks = [0.01, 0.1, 1]
)

clims = (0, 2.5)
clevels = range(clims...; step = 0.1)
cmap = :viridis

cb = Colorbar(pa[1,2],
    colormap = cgrad(cmap, length(clevels)-1; categorical=true),
    colorrange = clims,
    ticks = 0:4,
    ticksvisible = false,
    tickformat = values -> [rich("10", superscript("$(Int(z))")) for z in values],
    label = Ic_str
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

CL = @. leftboundary(R.*1u"μm", 1u"nM", Ū*1u"μm/s", T̄*1u"ms", D̄*1u"μm^2/s", a) |> u"μM"
lines!(ax_a, ustrip.(R), ustrip.(CL);
    linewidth=5, linestyle=:dash, color=:orange2, alpha=0.5,
)
CR = @. rightboundary(R.*1u"μm", 1u"nM", Ū*1u"μm/s", T̄*1u"ms", D̄*1u"μm^2/s", a) |> u"μM"
lines!(ax_a, ustrip.(R), ustrip.(CR);
    linewidth=5, linestyle=:dash, color=:orange2, alpha=0.5,
)

scatter!(ax_a,
    Rx, repeat([1*Cmin], length(Rx)),
    marker = :utriangle,
    markersize = 45,
    strokewidth = 1,
    strokecolor = :white,
    color = palette(:thermal, length(Rx))
)

scatter!(ax_a,
    repeat([1*Rmin], length(Cx)), Cx,
    marker = :rtriangle,
    markersize = 45,
    strokewidth = 1,
    strokecolor = :white,
    color = palette(:batlow, length(Cx))
)

xlims!(ax_a, (Rmin, Rmax))
ylims!(ax_a, (Cmin, Cmax))

##
fig
