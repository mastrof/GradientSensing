using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using CairoMakie
using PublicationFiguresMakie; set_theme!(Publication)

##
function plot_phycosphere_hein()
    datasets = collect_results(datadir("Hein"), rinclude=[r"phycosphere"], rexclude=[r"old"])
    for df in eachrow(datasets)
        f = jldopen(datadir("Hein", "RC.jld2"))
        R, Cₛ = f["R"], f["Cₛ"]
        S = df.S
        x = ustrip.(R)
        y = ustrip.(Cₛ)
        z = ustrip.(log10.(S ./ R))
        
        fig = Figure()
        cbar = Colorbar(fig[1,2], colormap=:viridis, colorrange=(0,2), label="log(S/R)")
        ax, plt = contourf(fig[1,1], x, y, z,
            levels = 0:0.25:2,
            xlabel = "R (μm)",
            ylabel = "Cₛ (μm)",
            xscale = log10,
            yscale = log10
        )

    #==
        contourf(R, Cₛ, log10.(S./R)', scale=:log10)
        plot!(
            xlab="R (μm)", ylab="Cₛ (μM)", cbartitle="log(S/R)",
            yticks=exp10.(-10:2:4)
        )
    ==#
        # parse dataset parameters
        #params = parse_savename(df.path)[2]
        # save both svg and png
        #savefig(plotsdir("Hein", savename("phycosphere", params, "svg")))
        # savefig(plotsdir("Hein", savename("phycosphere", params, "png")))
        #save(plotsdir("Hein", savename("phycosphere", params, "svg")))
    end
end

plot_phycosphere_hein()

##
function plot_phycosphere_hein_byDc()
    f = jldopen(datadir("Hein", "RC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]
    close(f)

    alldatasets = unpack_dataframe(
        collect_results(datadir("Hein"), rinclude=[r"phycosphere"], rexclude=[r"old"])
    )
    datasets = subset(alldatasets,
        :U => U -> U .== 50,
        :T => T -> T .== 100,
        :Π => Π -> Π .== 6,
    )
    sort!(datasets, [:Dc])

    plot(
        layout=(1,3), size=(1600,600),
        xlab="R (μm)", ylab="Cₛ (μM)", cbartitle="log(S/R)",
    )
    for (i,df) in enumerate(eachrow(datasets))
        S = df.S
        contourf!(subplot=i,
            ustrip.(R), ustrip.(Cₛ), log10.(S ./ R)',
            clims=(0,2), scale=:log10
        )
        if i ≠ 3
            plot!(subplot=i, colorbar=false)
        end
    end
    savefig(plotsdir("Hein", "phycosphere_varDc_T=100_U=50_Π=6.svg"))
end

plot_phycosphere_hein_byDc()
