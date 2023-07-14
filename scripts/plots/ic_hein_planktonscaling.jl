##
using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using CairoMakie
using PublicationFiguresMakie
set_theme!(Publication,
    Axis=(
        xticksvisible=false, yticksvisible=false,
        xminorticksvisible=false, yminorticksvisible=false,
    )
)

##
f = jldopen(datadir("HeinMod", "RC.jld2"), "r")
R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
close(f)

##
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

##
PER = [0.02, 0.1, 0.5]
r = R[0.5 .< R .< 60]
L = leakage_rate.(r .* 1u"μm", PER')
Cp = C.(r .* 1u"μm", L, 1u"nM", D̄*1u"μm^2/s")
Cp = ustrip.(Cp)

##
clims = (0, 3.0)
clevels = range(clims..., step=0.25)
cmap = :viridis

fig = Figure(resolution = SinglePlot());
Colorbar(fig[1,2],
    colormap = cmap,
    colorrange = clims,
    ticks = 0:4,
    label = "log(IC)"
)

ic = df.ic
ax, plt = contourf(fig[1,1], R, Cₛ, log10.(ic),
    colormap = cmap,
    levels = clevels,
)

ax.xscale = log10
ax.yscale = log10
ax.xlabel = "R (μm)"
ax.ylabel = "Cₛ (μM)"
ylims!(Cₛ[1], Cₛ[end])
xlims!(r[1], r[end])

band!(ax, r, Cp[:,1], Cp[:,3], color=RGBAf(0.96,0.9,0.9,0.35))
lines!(ax, r, Cp[:,2], linewidth=8, color=RGBAf(0.95,0.95,0.95,0.95))
fig
