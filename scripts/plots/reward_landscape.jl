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

##
R̄ = [0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0]
js = map(r -> findfirst(R .> r), R̄)

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
#save(plotsdir("HeinMod", "reward-vs-C.svg"), fig)

##
C̄ = [0.01, 0.015, 0.02, 0.04, 0.1]
js = map(c -> findfirst(Cₛ .> c), C̄)

C̄ = Cₛ[js]
Ī = ic[:,js]'
cmap = palette(:thermal, length(C̄))
fig = Figure(resolution = SinglePlot());
labels = string.(round.(C̄ .* 1e3, sigdigits=2))
ax, _ = series(fig[1,1], R, Ī, color=cmap, linewidth=9, labels=labels)
tightlimits!(ax)
ax.yscale = log10
ax.xscale = log10
ax.xlabel = "R (μm)"
ax.ylabel = "IC"
xlims!(ax, 0.5, 50)
leg = axislegend("Cₛ (nM)"; position=:rt)
#save(plotsdir("HeinMod", "reward-vs-R.svg"), fig)
