using DrWatson
@quickactivate :GradientSensing
using DataFrames
using JLD2
using Roots
using GLMakie
using PublicationFiguresMakie
arial_italic = "/usr/share/fonts/TTF/ariali.ttf"
set_theme!(Publication,
    Axis = (
        xminorticksvisible=false, yminorticksvisible=false,
        titlesize=32, titlefont=:bold, titlealign=:left,
    ),
    fonts = (
        regular="Arial",
        bold="Arial Bold",
        italic=arial_italic,
    ),
)

## load two sample molecular adsorption times
f = jldopen(datadir("Poisson", "larger.jld2"))
waitingtimes = f["waitingtimes"][10][50]
close(f)
adsorption_times = cumsum(waitingtimes) # sensory timescale T = 100ms
# split around midpoint
t1 = adsorption_times[adsorption_times .< 50u"ms"]
t2 = adsorption_times[adsorption_times .>= 50u"ms"]
t2 .-= t1[end]

## load the results of a KS simulation
R = 6.26
f = jldopen(datadir("Poisson", "KolmogorovSmirnov", savename("sensing", Dict(
    "C₀" => 1, # nM
    "Cₛ" => 0.25, # μM
    "Dc" => 500, # μm^2/s
    "N" => 500,
    "R" => R, # μm
    "T" => 100, # ms
    "U" => 50, # μm/s
    "Δt" => 0.0001, # ms
), "jld2")))
r = f["r"]
ksavg = f["ksavg"]
ks = f["ks"]
close(f)
ravg = mean(hcat(r...); dims=2)[:,1]
S = ustrip(ravg[findfirst(ksavg .>= 0.99)])


## figure
fig = Figure(size=TwoColumns(1.25))
# dummy space
Label(fig[1:2,2], "")
# ks ensemble
ax3 = Axis(fig[1:2,3];
    xlabel="distance from source (μm)", ylabel="detected a gradient",
    yticks=[0,1], yticksmirrored=false,
    ytickformat=(values -> [v == 0 ? "no" : "yes" for v in values]),
    ylabelcolor=cgrad(:Dark2)[4], yticklabelcolor=cgrad(:Dark2)[4],
    title="C"
)
ax4 = Axis(fig[1:2,3];
    ylabel="consensus fraction", yaxisposition=:right, yticks=[0,0.5,1],
    yticksmirrored=false,
    ylabelcolor=cgrad(:Dark2)[3], yticklabelcolor=cgrad(:Dark2)[3]
)
ylims!(ax3, -0.2, 1.2)
ylims!(ax4, -0.2, 1.2)
linkxaxes!(ax3, ax4)
hidespines!(ax4)
hidexdecorations!(ax4)
n = 500
for i in 1:n
    x = ustrip(r[i])
    y = @. ks[i] + (i-n/2)*0.0005
    scatter!(ax3, x, y; markersize=12, alpha=0.06, color=cgrad(:Dark2)[4])
end
scatterlines!(ax4, ustrip(ravg), ksavg;
    color=cgrad(:Dark2)[3], linewidth=6, markersize=24
)
errorbars!(ax4, ustrip(ravg), ksavg,
    [std([ks[j][i] for j in eachindex(ks)])/sqrt(length(ksavg)) for i in eachindex(ksavg)];
    color=cgrad(:Dark2)[3], linewidth=6
)
vlines!(ax3, [R]; linestyle=:dash, color=:gray, linewidth=2)
vlines!(ax3, [S];
    linestyle=:dash, color=cgrad(:Dark2)[2], linewidth=4
)
# adsorption times
ax1 = Axis(fig[1,1]; xlabel="time (ms)", title="A")
hideydecorations!(ax1)
vlines!(ax1, ustrip.(t1); color=cgrad(:Dark2)[1])
vlines!(ax1, ustrip.(t2 .+ t1[end]); color=cgrad(:Dark2)[2])
vlines!(ax1, [50]; linewidth=5, linestyle=:dash, color=:black)
# ecdf
ax2 = Axis(fig[2,1];
    xlabel="waiting time (ms)", ylabel="cumulative distribution", title="B"
)
ecdfplot!(ax2, diff(ustrip.(t1)); linewidth=4)
ecdfplot!(ax2, diff(ustrip.(t2)); linewidth=4)
fig
