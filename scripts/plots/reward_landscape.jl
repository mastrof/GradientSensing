##
using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using CairoMakie
using PublicationFiguresMakie
set_theme!(Publication,
    Axis = (
        xminorticksvisible = false,
        yminorticksvisible = false,
    )
)

##
alldatasets = unpack_dataframe(
    collect_results(datadir("HeinMod"); rinclude=[r"IC"])
)

##
D̄ = 500
Ū = 50
T̄ = 100
Π̄ = 6
ic = subset(alldatasets,
    :Dc => Dc -> Dc .== D̄,
    :U => U -> U .== Ū,
    :T => T -> T .== T̄,
    :Π => Π -> Π .== Π̄
)[1,:].ic

##
f = jldopen(datadir("HeinMod", "RC.jld2"), "r")
R = ustrip.(f["R"])
Cₛ = ustrip.(f["Cₛ"])
close(f)
R̄ = [0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0]
js = map(r -> findfirst(R .> r), R̄)

##
R̄ = R[js]
Ī = ic[js,:]
cmap = palette(:thermal, length(R̄))
fig = Figure(resolution = SinglePlot());
labels = string.(round.(R̄, sigdigits=2))
ax, _ = series(fig[1,1], Cₛ, Ī, color=cmap, linewidth=9, labels=labels)
tightlimits!(ax)
ax.yscale = log10
ax.xscale = log10
ax.xlabel = "Cₛ (μM)"
ax.ylabel = "IC"
xlims!(ax, 3e-3, Cₛ[end])
leg = axislegend("R (μm)"; position=:lt)
save(plotsdir("HeinMod", "reward-vs-C.svg"), fig)
