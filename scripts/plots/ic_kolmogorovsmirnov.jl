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
        datadir("newPoisson","KolmogorovSmirnov"),
        rinclude=[r"IC"]
    )
    for df in eachrow(datasets)
        f = jldopen(datadir("newPoisson", "RC.jld2"))
        R, Cₛ = f["R"], f["Cₛ"]
        R = R[1:end-1]
        ic = df.ic

        heatmap(R, Cₛ, log10.(ic)', scale=:log10)
        plot!(xlab="R (μm)", ylab="Cₛ (μM)", cbartitle="log(IC)")

        # parse dataset parameters
        params = parse_savename(df.path)[2]
        # save both svg and png
        savefig(plotsdir("newPoisson", "KolmogorovSmirnov", savename("IC", params, "svg")))
        savefig(plotsdir("newPoisson", "KolmogorovSmirnov", savename("IC", params, "png")))
    end
end

plot_ic_ks()
