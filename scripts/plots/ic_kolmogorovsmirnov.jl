using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using Plots; Plots.pyplot()
Plots.default(
    thickness_scaling=1.5, frame=:border, grid=false,
    bgcolorlegend=:transparent, fgcolorlegend=:transparent,
    palette=:Dark2
)

function plot_ic_ks()
    datasets = collect_results(
        datadir("PoissonSampling","KolmogorovSmirnov"),
        rinclude=[r"IC"]
    )
    for df in eachrow(datasets)
        # load R and L
        f = jldopen(datadir("PoissonSampling", "paramspaceRC.jld2"))
        R, Cₛ = f["R"], f["Cₛ"]
        ic = df.ic
        
        heatmap(R, Cₛ, log10.(ic)', scale=:log10)
        plot!(
            xlab="R (μm)", ylab="Cₛ (μM)", cbartitle="log(IC)",
            yticks=exp10.(-4:0)
        )

        # parse dataset parameters
        params = parse_savename(df.path)[2]
        # save both svg and png
        savefig(plotsdir("PoissonSampling", "KolmogorovSmirnov", savename("IC", params, "svg")))
        savefig(plotsdir("PoissonSampling", "KolmogorovSmirnov", savename("IC", params, "png")))
    end
end

plot_ic_ks()