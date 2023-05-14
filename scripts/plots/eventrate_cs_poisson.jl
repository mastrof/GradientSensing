using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using Plots; Plots.pyplot()
Plots.default(
    thickness_scaling=1.5, frame=:border, grid=false,
    bgcolorlegend=:transparent, fgcolorlegend=:transparent,
    palette=:Dark2
)

function plot_eventrate_poisson()
    f = jldopen(datadir("Poisson", "RC.jld2"))
    Cₛ = f["Cₛ"]
    E = eventrate.(Cₛ)
    plot(Cₛ, E, scale=:log10, lw=5, leg=false)
    plot!(xlab="Cₛ (μM)", ylab="E (1/s)")
    savefig(plotsdir("Poisson", "eventrate-vs-Cs.svg"))
end

plot_eventrate_poisson()
