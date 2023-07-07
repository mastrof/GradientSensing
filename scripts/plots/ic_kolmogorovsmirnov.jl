##
using DrWatson
@quickactivate :GradientSensing
using JLD2, DataFrames
using CairoMakie
using PublicationFiguresMakie
set_theme!(Publication,
    Axis=(
        xticksvisible=false, yticksvisible=false,
        xminorticksvisible=false, yminorticksvisible=false,
    )
)


##
function makeplot(datasets, R, Cₛ)
    clims = (0, 4.5)
    clevels = range(clims..., step=0.25)
    cmap = :viridis

    fig = Figure(resolution = TwoColumns(1));
    Colorbar(fig[1,4],
        colormap = cmap,
        colorrange = clims,
        ticks = 0:4,
        label = "log(IC)"
    )
    labels = ["A", "B", "C"]

    for (i,df) in enumerate(eachrow(datasets))
        ic = df.ic
        ax, _ = heatmap(fig[1,i], R, Cₛ, log10.(ic),
            colormap = cmap,
            colorrange = clims,
        )
        tightlimits!(ax)
        ax.xscale = log10
        ax.yscale = log10
        ax.xlabel = "R (μm)"
        if i == 1
            ax.ylabel = "Cₛ (μM)"
        end
        Label(fig[1, i, TopLeft()], labels[i], halign=:right)
    end
    fig
end

##
#==
function plot_phycosphere_hein_byDc()
    f = jldopen(datadir("Poisson", "RC.jld2"))
    R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
    close(f)

    alldatasets = unpack_dataframe(
        collect_results(datadir("Poisson", "KolmogorovSmirnov"), rinclude=[r"IC"])
    )
    Ū = 50
    T̄ = 100
    datasets = subset(alldatasets,
        :U => U -> U .== Ū,
        :T => T -> T .== T̄,
    )
    sort!(datasets, [:Dc])
    
    fig = makeplot(datasets, R, Cₛ);
    for (i,df) in enumerate(eachrow(datasets))
        ax = contents(fig[1,i])[1]
        ax.title = "Dc = $(df.Dc)μm²/s"
    end

    save(plotsdir("Poisson", "KolmogorovSmirnov", "IC_varDc_T=$(T̄)_U=$(Ū).svg"), fig)
end

plot_phycosphere_hein_byDc()
==#

##
#==
function plot_phycosphere_hein_byU()
    f = jldopen(datadir("Poisson", "RC.jld2"))
    R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
    close(f)

    alldatasets = unpack_dataframe(
        collect_results(datadir("Poisson", "KolmogorovSmirnov"), rinclude=[r"IC"])
    )
    D̄ = 500
    T̄ = 100
    datasets = subset(alldatasets,
        :Dc => D -> D .== D̄,
        :T => T -> T .== T̄,
    )
    sort!(datasets, [:U])
    
    fig = makeplot(datasets, R, Cₛ)
    for (i,df) in enumerate(eachrow(datasets))
        ax = contents(fig[1,i])[1]
        ax.title = "U = $(df.U)μm/s"
    end

    save(plotsdir("Poisson", "KolmogorovSmirnov", "IC_varU_Dc=$(D̄)_T=$(T̄).svg"), fig)
end

plot_phycosphere_hein_byU()
==#

##
function plot_phycosphere_hein_byT()
    f = jldopen(datadir("Poisson", "RC.jld2"))
    R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
    close(f)

    alldatasets = unpack_dataframe(
        collect_results(datadir("Poisson", "KolmogorovSmirnov"), rinclude=[r"IC"])
    )
    D̄ = 500
    Ū = 50
    datasets = subset(alldatasets,
        :Dc => D -> D .== D̄,
        :U => U -> U .== Ū,
    )
    sort!(datasets, [:T])
    
    fig = makeplot(datasets, R, Cₛ);
    for (i,df) in enumerate(eachrow(datasets))
        ax = contents(fig[1,i])[1]
        ax.title = "T = $(df.T)ms"
    end

    save(plotsdir("Poisson", "KolmogorovSmirnov", "IC_varT_Dc=$(D̄)_U=$(Ū).svg"), fig)
end

plot_phycosphere_hein_byT()
