using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using Plots; Plots.pyplot()
Plots.default(
    thickness_scaling=1.5, frame=:border, grid=false,
    bgcolorlegend=:transparent, fgcolorlegend=:transparent,
    palette=:Dark2
)

function snr_kolmogorovsmirnov()
    datasets = collect_results(
        datadir("Poisson", "KolmogorovSmirnov"),
        rinclude = [r"SNR"]
    )
    for df in eachrow(datasets)
        f = jldopen(datadir("Poisson", "RC.jld2"))
        @unpack R, Cₛ = f
        SNR = df.SNR

        heatmap(R, Cₛ, SNR', scale=:log10)
        contour!(R, Cₛ, SNR', levels=[0.25, 0.5, 1], lw=5, c=[:green2, :cyan, :pink])
        plot!(xlab="R (μm)", ylab="Cₛ (μM)", cbartitle="SNR @ phycosphere boundary")

        params = parse_savename(df.path)[2]
        savefig(plotsdir("Poisson", "KolmogorovSmirnov", savename("snr", params, "svg")))
        savefig(plotsdir("Poisson", "KolmogorovSmirnov", savename("snr", params, "png")))
    end
end

snr_kolmogorovsmirnov()
