using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using Plots; Plots.pyplot()
Plots.default(
    thickness_scaling=1.5, frame=:border, grid=false,
    bgcolorlegend=:transparent, fgcolorlegend=:transparent,
    palette=:Dark2
)

function plot_eventrate_hein()
    datasets = collect_results(datadir("Hein"), rinclude=[r"phycosphere"])
    for df in eachrow(datasets)
        f = jldopen(datadir("Hein", "RC.jld2"))
        R, Cₛ = f["R"], f["Cₛ"]
        S = df.S
        C_at_S = C.(S, R', Cₛ, C₀)
        E = eventrate.(C_at_S)
        E₀ = eventrate(C₀)

        contourf(R, Cₛ, log10.(E'./E₀), scale=:log10)
        plot!(xlab="R (μm)", ylab="Cₛ (μM)", cbartitle="log(E/E₀)")

        params = parse_savename(df.path)[2]
        savefig(plotsdir("Hein", savename("eventrate-at-S", params, "svg")))
    end
end

plot_eventrate_hein()
