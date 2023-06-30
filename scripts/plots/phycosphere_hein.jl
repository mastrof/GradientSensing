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
        titlefont=:regular, titlesize=24,
    )
)

##
function plot_phycosphere_hein(df)
    f = jldopen(datadir("Hein", "RC.jld2"))
    R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
    S = ustrip.(df.S)
    
    fig = Figure();
    cbar = Colorbar(fig[1,2], colormap=:viridis, colorrange=(0,2.5), label="log(S/R)")
    ax, plt = contourf(fig[1,1], R, Cₛ, log10.(S./R), levels = 0:0.25:2.5)

    xlims!(ax, low=R[1], high=R[end])
    ylims!(ax, low=Cₛ[1], high=Cₛ[end])
    ax.xscale = log10
    ax.yscale = log10
    ax.xlabel = "R (μm)"
    ax.ylabel = "Cₛ (μM)"
    ax.title = gettitle(df.path)
    fig
end

function gettitle(fname)
    config = parse_savename(fname)[2]
    title = join(
        [
            "Dc = $(config["Dc"])μm²/s",
            "U = $(config["U"])μm/s",
            "T = $(config["T"])ms",
            "Π = $(config["Π"])"
        ]
        , ", "
    )
    return title
end

##
function plot_phycosphere_hein()
    datasets = collect_results(datadir("Hein"), rinclude=[r"phycosphere"], rexclude=[r"old"])
    for df in eachrow(datasets)
        fig = plot_phycosphere_hein(df)
        config = parse_savename(df.path)[2]
        save(plotsdir("Hein", savename("phycosphere", config, "svg")), fig)
    end
end

plot_phycosphere_hein()

##
function plot_phycosphere_hein_byDc()
    f = jldopen(datadir("Hein", "RC.jld2"))
    R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
    close(f)

    alldatasets = unpack_dataframe(
        collect_results(datadir("Hein"), rinclude=[r"phycosphere"], rexclude=[r"old"])
    )
    Ū = 50
    T̄ = 100
    Π̄ = 6
    datasets = subset(alldatasets,
        :U => U -> U .== Ū,
        :T => T -> T .== T̄,
        :Π => Π -> Π .== Π̄,
    )
    sort!(datasets, [:Dc])
    
    fig = Figure(resolution = TwoColumns(1));
    cbar = Colorbar(fig[1,4], colormap = :viridis, colorrange = (0, 2.5), label = "log(S/R)")
    labels = ["A", "B", "C"]
    for (i,df) in enumerate(eachrow(datasets))
        S = ustrip.(df.S)
        ax, plt = contourf(fig[1,i], R, Cₛ, log10.(S./R), colormap = :viridis, levels = 0:0.25:2.5)
        tightlimits!(ax)
        ax.xscale = log10
        ax.yscale = log10
        ax.xlabel = "R (μm)"
        if i == 1
            ax.ylabel = "Cₛ (μM)"
        end
        ax.title = "Dc = $(df.Dc)μm²/s"
        Label(fig[1,i,TopLeft()], labels[i],
            fontsize = 32,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end

    save(plotsdir("Hein", "phycosphere_varDc_T=$(T̄)_U=$(Ū)_Π=$(Π̄).svg"), fig)
end

plot_phycosphere_hein_byDc()


##
function plot_phycosphere_hein_byU()
    f = jldopen(datadir("Hein", "RC.jld2"))
    R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
    close(f)

    alldatasets = unpack_dataframe(
        collect_results(datadir("Hein"), rinclude=[r"phycosphere"], rexclude=[r"old"])
    )
    D̄ = 500
    T̄ = 100
    Π̄ = 6
    datasets = subset(alldatasets,
        :Dc => D -> D .== D̄,
        :T => T -> T .== T̄,
        :Π => Π -> Π .== Π̄,
    )
    sort!(datasets, [:U])
    
    fig = Figure(resolution = TwoColumns(1));
    cbar = Colorbar(fig[1,4], colormap = :viridis, colorrange = (0, 2.5), label = "log(S/R)")
    labels = ["A", "B", "C"]
    for (i,df) in enumerate(eachrow(datasets))
        S = ustrip.(df.S)
        ax, plt = contourf(fig[1,i], R, Cₛ, log10.(S./R), colormap = :viridis, levels = 0:0.25:2.5)
        tightlimits!(ax)
        ax.xscale = log10
        ax.yscale = log10
        ax.xlabel = "R (μm)"
        if i == 1
            ax.ylabel = "Cₛ (μM)"
        end
        ax.title = "U = $(df.U)μm/s"
        Label(fig[1,i,TopLeft()], labels[i],
            fontsize = 32,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end

    save(plotsdir("Hein", "phycosphere_varU_Dc=$(D̄)_T=$(T̄)_Π=$(Π̄).svg"), fig)
end

plot_phycosphere_hein_byU()

##
function plot_phycosphere_hein_byT()
    f = jldopen(datadir("Hein", "RC.jld2"))
    R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
    close(f)

    alldatasets = unpack_dataframe(
        collect_results(datadir("Hein"), rinclude=[r"phycosphere"], rexclude=[r"old"])
    )
    D̄ = 500
    Ū = 50
    Π̄ = 6
    datasets = subset(alldatasets,
        :Dc => D -> D .== D̄,
        :U => U -> U .== Ū,
        :Π => Π -> Π .== Π̄,
    )
    sort!(datasets, [:T])
    
    fig = Figure(resolution = TwoColumns(1));
    cbar = Colorbar(fig[1,4], colormap = :viridis, colorrange = (0, 2.5), label = "log(S/R)")
    labels = ["A", "B", "C"]
    for (i,df) in enumerate(eachrow(datasets))
        S = ustrip.(df.S)
        ax, plt = contourf(fig[1,i], R, Cₛ, log10.(S./R), colormap = :viridis, levels = 0:0.25:2.5)
        tightlimits!(ax)
        ax.xscale = log10
        ax.yscale = log10
        ax.xlabel = "R (μm)"
        if i == 1
            ax.ylabel = "Cₛ (μM)"
        end
        ax.title = "T = $(df.T)ms"
        Label(fig[1,i,TopLeft()], labels[i],
            fontsize = 32,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end

    save(plotsdir("Hein", "phycosphere_varT_Dc=$(D̄)_U=$(Ū)_Π=$(Π̄).svg"), fig)
end

plot_phycosphere_hein_byT()
