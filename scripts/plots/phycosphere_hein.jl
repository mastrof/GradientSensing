using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using Plots; Plots.pyplot()
Plots.default(
    thickness_scaling=1.5, frame=:border, grid=false,
    bgcolorlegend=:transparent, fgcolorlegend=:transparent,
    palette=:Dark2
)

function plot_phycosphere_hein()
    datasets = collect_results(datadir("Hein"), rinclude=[r"phycosphere"])
    for df in eachrow(datasets)
        f = jldopen(datadir("Hein", "RC.jld2"))
        R, Cₛ = f["R"], f["Cₛ"]
        S = df.S

        contourf(R, Cₛ, log10.(S./R)', scale=:log10)
        plot!(
            xlab="R (μm)", ylab="Cₛ (μM)", cbartitle="log(S/R)",
            yticks=exp10.(-10:2:4)
        )

        # parse dataset parameters
        params = parse_savename(df.path)[2]
        # save both svg and png
        savefig(plotsdir("Hein", savename("phycosphere", params, "svg")))
        # savefig(plotsdir("Hein", savename("phycosphere", params, "png")))
    end
end

plot_phycosphere_hein()
