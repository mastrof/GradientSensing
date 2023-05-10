using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using Plots; Plots.pyplot()
Plots.default(
    thickness_scaling=1.5, frame=:border, grid=false,
    bgcolorlegend=:transparent, fgcolorlegend=:transparent,
    palette=:Dark2
)

function plot_phycosphere_ks()
    datasets = collect_results(
        datadir("PoissonSampling", "KolmogorovSmirnov"),
        rinclude=[r"phycosphere"]
    )
    for df in eachrow(datasets)
        # load R and Cₛ
        f = jldopen(datadir("PoissonSampling", "paramspaceRC.jld2"))
        R, Cₛ = f["R"], f["Cₛ"]
        S = df.S

        heatmap(R, Cₛ, log10.(S./R)', scale=:log10)
        plot!(
            xlab="R (μm)", ylab="Cₛ (pmol/s)", cbartitle="log(S/R)",
            yticks=exp10.(-9:2:5)
        )

        # parse dataset parameters
        params = parse_savename(df.path)[2]
        # save both svg and png
        savefig(plotsdir("PoissonSampling", "KolmogorovSmirnov", savename("phycosphere", params, "svg")))
        savefig(plotsdir("PoissonSampling", "KolmogorovSmirnov", savename("phycosphere", params, "png")))
    end
end

plot_phycosphere_ks()
