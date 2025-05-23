##
using DrWatson
@quickactivate :GradientSensing
using DataFrames
using JLD2
using StatsBase
using Interpolations
using Roots
using Unitful
using GLMakie
using ColorSchemes
palette(colors::Symbol, n::Int) = get(colorschemes[colors], range(0,1;length=n))

##
function community_structure(α, Ntot, nclasses; r_min=0.5, r_max=70)
    V_min, V_max = @. 4π / 3 * (r_min^3, r_max^3) # μm³
    V_edges = exp2.(range(log2(V_min), log2(V_max); length=nclasses + 1))
    V = midpoints(V_edges)
    r_edges = @. (3 * V_edges / (4π))^(1 / 3)
    r = @. (3 * V / (4π))^(1 / 3) # μm
    N = @. V^α
    N = N ./ sum(N) .* Ntot # cells / mL
    return r, r_edges, N
end

function encounter_times!(df, placeholder, N, λ, U, T, PER; Dc=500, Π=6, C0=1u"nM")
    Te = @. time_to_encounter(r+ustrip(a), λ, U, N) / 3600 # hours
    r_units = r .* 1u"μm"
    re_units = r_edges .* 1u"μm"
    L = leakage_rate.(re_units, PER)
    Dc_units = Dc * 1u"μm^2/s"
    U_units = U * 1u"μm/s"
    Cp_units = C.(re_units, L, 1u"nM", Dc_units)
    Se = @. HeinModRadius(re_units, Cp_units, C0, T * 1u"ms", Dc_units, U_units, Π) |> u"μm" |> ustrip
    #S = @. HeinModRadius(r_units, Cp_units, C0, T * 1u"ms", Dc_units, U_units, Π) |> u"μm" |> ustrip
    S = [sqrt(Se[i-1]*Se[i]) for i in eachindex(Se)[2:end]]
    Ice = @. IC(λ, r_edges + ustrip(a), Se + ustrip(a))
    Ic = [sqrt(Ice[i-1]*Ice[i]) for i in eachindex(Ice)[2:end]]
    #Ic = @. IC(λ, r + ustrip(a), S + ustrip(a))
    newdf = DataFrame(
        U=U,
        T=T,
        R=[r],
        S=[S],
        Ic=[Ic],
        PER=PER,
        Te=[Te ./ Ic],
    )
    append!(df, newdf)
end

function dfsetup()
    DataFrame(
        U=Float64[],
        T=Float64[],
        R=Vector{Float64}[],
        S=Vector{Float64}[],
        Ic=Vector{Float64}[],
        PER=Float64[],
        Te=Vector{Float64}[]
    )
end

function time_to_encounter(R, λ, U, N)
    Γ = encounter_kernel(R, λ, U) # μm³/s
    Γ *= 1e-12 # mL / s
    T = 1 / (Γ * N) # s
    return T
end

function encounter_kernel(R, λ, U)
    kn = 4λ / 3R
    β = (kn / (1 + kn))
    return β * π * U * R^2
end

function find_r_prc(N, r, p)
    Ntot = sum(N)
    f = cumsum(N) ./ Ntot
    fitp = linear_interpolation(r, f)
    rp = find_zero(x -> fitp(x) - p, r[2])
    return rp
end

## Parameter values used throughout the entire code
Ntot_oligotrophic = 1e5 # total phytoplankton abundance (cell/mL)
Ntot_productive = 2e5 # (cell/mL)
n = 17 # number of size classes in community size spectrum
r_min, r_max = 0.5, 70.0 # extremal values of phytoplankton size
V_min, V_max = @. 4π / 3 * (r_min^3, r_max^3) # μm³
V_edges = exp2.(range(log2(V_min), log2(V_max); length=n+1))
V = midpoints(V_edges)
r_edges = @. (3 * V_edges / (4π))^(1 / 3)
r = @. (3 * V / (4π))^(1 / 3) # μm
Us = [25 60] # swimming speeds
Ts = [250 100] # sensory timescales
λ = 30 # bacterial correlation length (μm)
PER = [0.02, 0.1, 0.4] # values of percent extracellular release

##

# define a range of volumes uniformly distributed in log2 scale
α_oligotrophic = -1
α_productive = -0.75
r, r_edges, N_oligotrophic = community_structure(α_oligotrophic, Ntot_oligotrophic, n; r_min, r_max)
_, _, N_productive = community_structure(α_productive, Ntot_productive, n; r_min, r_max)
N = Dict(
    :oligotrophic => N_oligotrophic,
    :productive => N_productive
)
prc = Dict(
    :oligotrophic => [find_r_prc(N[:oligotrophic], r, p) for p in [0.95, 0.99]],
    :productive => [find_r_prc(N[:productive], r, p) for p in [0.95, 0.99]]
)

## evaluate times to encounters
df = Dict(
    :oligotrophic => dfsetup(),
    :productive => dfsetup()
)
for environment in [:oligotrophic, :productive], per in PER, (T, U) in zip(Ts, Us)
    C0 = environment == :oligotrophic ? 1u"nM" : 100u"nM"
    encounter_times!(df[environment], r, N[environment], λ, U, T, per; C0)
end


## Initialize figure layout
fig = Figure(resolution=(1600,720).*1)
panela = fig[1:2, 1] = GridLayout()
panelb = fig[3:4, 1] = GridLayout()
panelc = fig[1:3, 2] = GridLayout()
paneld = fig[1:3, 3] = GridLayout()
cdlegend = fig[4, 2:3]

getrgb(c) = (c.r, c.g, c.b)
alphaize(c, alpha) = RGBAf(getrgb(c)..., alpha)

# panel C -- search times in oligotrophic environment
gdf = Dict(k => groupby(df[k], :T) for k in keys(df))

ax_oligo = Axis(panelc[1, 1];
    xlabel=rich(rich("R"; font=:italic), " (μm)"),
    ylabel=rich(rich("T", subscript("e"); font=:italic), " (hours)"),
    yscale=log10,
    xscale=log10,
    yticks=exp10.(0:3),
    xticks=[1, 3, 9, 27],
    title="Oligotrophic waters"
)
xlims!(ax_oligo, (0.45, 75))
ylims!(ax_oligo, (5.5, 2e3))
colors = palette(:Dark2_8, 8)[[1, 3]]

bandalpha = 0.35
curvelab(U, T) = rich(
    rich("U"; font=:italic), "=$(Int(U))μm/s, ",
    rich("T"; font=:italic), "=$(Int(T))ms",
)
labs = @. rich(curvelab(Us, Ts); color=colors)

for (i, g) in enumerate(gdf[:oligotrophic])
    t_low = g[g.PER.==PER[1], :Te][1]
    t_mid = g[g.PER.==PER[2], :Te][1]
    t_hi = g[g.PER.==PER[3], :Te][1]
    craw = colors[i]
    c = alphaize(craw, bandalpha)
    band!(ax_oligo, r, t_low, t_hi; color=c)
    scatterlines!(ax_oligo, r, t_mid; color=craw,
        linewidth=8,
        markersize=24,
        marker=i == 1 ? :circle : :rect
        #label=labs[:, i]
    )
end

# random search times for reference
Te_rnd_oligo = @. time_to_encounter(r+ustrip(a), λ, Us, N[:oligotrophic]) / 3600
lines!(ax_oligo, r, Te_rnd_oligo[:, 1];
    linestyle=:dash,
    linewidth=3,
    color=alphaize(colors[1], 0.7)
)
lines!(ax_oligo, r, Te_rnd_oligo[:, 2];
    linestyle=:dash,
    linewidth=3,
    color=alphaize(colors[2], 0.7)
)

hlines!(ax_oligo, [24, 24 * 7, 24 * 30]; linestyle=:dash, color=:black, label=false)
vlines!(ax_oligo, prc[:oligotrophic]; linestyle=:dot, color=:black, label=false)
text!(ax_oligo,
    [25, 25, 25], [24, 24 * 7, 24 * 30];
    text=["1 day", "1 week", "1 month"],
    fontsize=22,
    align=(:left, :top)
)
text!(ax_oligo,
    prc[:oligotrophic] .* [0.55, 1.05], [1300];
    text=["95%", "99%"],
    fontsize=22,
    align=(:left, :center)
)


# panel D -- search times in productive environment
gdf = Dict(k => groupby(df[k], :T) for k in keys(df))

ax_prod = Axis(paneld[1, 1];
    xlabel=rich(rich("R"; font=:italic), " (μm)"),
    yscale=log10,
    xscale=log10,
    yticks=exp10.(0:3),
    yticklabelsvisible=false,
    xticks=[1, 3, 9, 27],
    title="Productive waters"
)
xlims!(ax_prod, (0.45, 75))
ylims!(ax_prod, (5.5, 2e3))
colors = palette(:Dark2_8, 8)[[1, 3]]

bandalpha = 0.35
curvelab(U, T) = rich(
    rich("U"; font=:italic), " = $(Int(U)) μm/s    ",
    rich("T"; font=:italic), " = $(Int(T)) ms",
)
labs = [rich(l; color=c) for (l, c) in zip(curvelab.(Us, Ts), colors)]

tactic = []
tacticbands = []
legtitles = []
for (i, g) in enumerate(gdf[:productive])
    t_low = g[g.PER.==PER[1], :Te][1]
    t_mid = g[g.PER.==PER[2], :Te][1]
    t_hi = g[g.PER.==PER[3], :Te][1]
    s = lines!(ax_prod, [r[1]], [t_mid[1]]; color=:white)
    push!(legtitles, s)
    craw = colors[i]
    c = alphaize(craw, bandalpha)
    s = band!(ax_prod, r, t_low, t_hi; color=c)
    push!(tacticbands, s)
    s = scatterlines!(ax_prod, r, t_mid; color=craw,
        linewidth=8,
        markersize=24,
        marker=i == 1 ? :circle : :rect
        #label=labs[:, i]
    )
    push!(tactic, s)
end

# random search times for reference
Te_rnd_prod = @. time_to_encounter(r+ustrip(a), λ, Us, N[:productive]) / 3600
random_1 = lines!(ax_prod, r, Te_rnd_prod[:, 1];
    linestyle=:dash,
    linewidth=3,
    color=alphaize(colors[1], 0.7)
)
random_2 = lines!(ax_prod, r, Te_rnd_prod[:, 2];
    linestyle=:dash,
    linewidth=3,
    color=alphaize(colors[2], 0.7)
)

hlines!(ax_prod, [24, 24 * 7, 24 * 30]; linestyle=:dash, color=:black, label=false)
vlines!(ax_prod, prc[:productive]; linestyle=:dot, color=:black, label=false)
#==
text!(ax_prod,
    [25, 0.5, 0.5], [24, 24 * 7, 24 * 30];
    text=["1 day", "1 week", "1 month"],
    fontsize=22,
    align=(:left, :bottom)
)
==#
text!(ax_prod,
    prc[:productive] .* [0.55, 1.05], [1300];
    text=["95%", "99%"],
    fontsize=22,
    align=(:left, :center)
)

# legend for panels C and D
legendlabels = [
    rich(labs[1]; fontsize=24),
    "Random encounters (PER = 0%)",
    "Chemotactic encounters (PER = 10%)",
    "Chemotactic encounters (2% < PER < 40%)"
]
legendlabels = repeat(legendlabels, 2)
legendlabels[5] = rich(labs[2]; fontsize=24)
Legend(cdlegend,
    [
        legtitles[1], random_1, tactic[1], tacticbands[1],
        legtitles[2], random_2, tactic[2], tacticbands[2]
    ],
    legendlabels,
    patchsize=(50, 25),
    nbanks=4,
    orientation=:horizontal,
    colgap=24,
    rowgap=6,
    valign=:bottom,
    margin=(0, 0, -10, 0)
)



# panel B -- community structures
ax_communities = Axis(panelb[1, 1];
    xlabel=rich(rich("R"; font=:italic), " (μm)"),
    ylabel=rich(rich("N"; font=:italic), " (cells/mL)"),
    xscale=log10,
    yscale=log10,
    xticks=[1, 3, 9, 27],
    yticks=exp10.(1:2:5),
    xtickformat=xs -> string.(Int.(xs)),
    ytickformat=ys -> [rich("10", superscript("$(Int(log10(y)))")) for y in ys],
    title="Community structure"
)
ylims!(ax_communities, (1, 1e5))

color_communities = palette(:tab10, 10)[1:2]

scatterlines!(ax_communities,
    r, N[:oligotrophic],
    linewidth=8,
    marker=:utriangle,
    markersize=32,
    color=color_communities[1],
    label=rich("Oligotrophic (10", superscript("5"), " cells/mL)")
)
scatterlines!(ax_communities,
    r, N[:productive],
    linewidth=8,
    marker=:dtriangle,
    markersize=32,
    color=color_communities[2],
    label=rich("Productive (2×10", superscript("5"), " cells/mL)")
)
axislegend(ax_communities; position=:rt, patchsize=(60, 20))


fig



## Load data for panel A
f = jldopen(datadir("HeinMod", "RC.jld2"), "r")
R, Cₛ = ustrip.(f["R"]), ustrip.(f["Cₛ"])
close(f)

Ū = 50
T̄ = 100
D̄ = 500
Π̄ = 6

alldatasets = unpack_dataframe(
    collect_results(datadir("HeinMod"), rinclude=[r"IC"])
)
df = subset(alldatasets,
    :U => U -> U .== Ū,
    :T => T -> T .== T̄,
    :Dc => Dc -> Dc .== D̄,
    :Π => Π -> Π .== Π̄
)[1, :]

## Phytoplankton scaling
PER = [0.02, 0.1, 0.5]
#r = R[0.5 .< R .< 60]
L = leakage_rate.(r .* 1u"μm", PER')
Cp = C.(r .* 1u"μm", L, 1u"nM", D̄ * 1u"μm^2/s")
Cp = ustrip.(Cp)


## Some global constants
Rmin, Rmax = r[1], r[end]
Cmin, Cmax = 1.5e-3, 1.0
cmap = :viridis
clims = (0, 2.5)
clevels = range(clims...; step=0.1)

Ic_str = rich("I", subscript("c"); font=:italic)

##
cb = Colorbar(panela[1, 2],
    colormap=cgrad(cmap, length(clevels) - 1; categorical=true),
    colorrange=clims,
    ticks=0:4,
    ticksvisible=false,
    tickformat=values -> [rich("10", superscript("$(Int(z))")) for z in values],
    label=Ic_str,
)

## IC landscape
ax = Axis(panela[1, 1],
    xscale=log10,
    yscale=log10,
    xticks=[1, 3, 9, 27],
    yticks=[0.01, 0.1, 1],
    xlabel=rich(rich("R"; font=:italic), " (μm)"),
    ylabel=rich(rich("C", subscript("s"); font=:italic), " (μM)"),
    title="Carbon-size scaling"
)
ylims!(Cmin, Cmax)
xlims!(Rmin, Rmax)
#colgap!(panela, 10)

ic = df.ic
contourf!(ax, R, Cₛ, log10.(ic),
    colormap=cmap,
    levels=clevels,
)

contour!(ax, R, Cₛ, log10.(df.ic),
    levels=[clevels[findfirst(clevels .> 0)]],
    linewidth=8,
    color=:white
)

band!(ax, r, Cp[:, 1], Cp[:, 3], color=RGBAf(1, 0.5, 0.7, 0.25))
lines!(ax, r, Cp[:, 2], linewidth=8, color=RGBAf(1, 0.5, 0.7, 0.95))


##
fig
