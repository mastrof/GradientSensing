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
        titlesize=32,
    ),
    fonts = (
        regular="Arial",
        bold="Arial Bold",
        italic=arial_italic,
    ),
)

## concentration field
Cdiff(r, R, C0, Cs) = (C0 + Cs*R / r) |> u"μM"
∇Cdiff(r, R, C0, Cs) = (Cs*R / r^2) |> u"μM/μm"
## sensing
snr(U, C, ∇C, r, R, C0, Cs, T, Dc, a) = signal(U,∇C,r,R,C0,Cs)/noise(U,C,r,R,C0,Cs,T,Dc,a)
signal(U, ∇C, r, R, C0, Cs) = U * ∇C(r,R,C0,Cs) |> u"μM/s"
noise(U, C, r, R, C0, Cs, T, Dc, a) = 6 * σ0(C,r,R,C0,Cs,T,Dc,a) |> u"μM/s"
σ0(C,r,R,C0,Cs,T,Dc,a) = sqrt(3*C(r,R,C0,Cs) / (π*a*Dc*T^3*Unitful.Na))
ξ(U, r, T) = 1 / (1 - exp(-(r/(U*T))^(3/2)))

## parameter space
f = jldopen(datadir("HeinMod", "RC.jld2"), "r")
radii, concentrations = f["R"], f["Cₛ"]
close(f)

U = 50u"μm/s"
T = 100u"ms"
Dc = 500u"μm^2/s"
a = 0.5u"μm"

S(C0) = map(Iterators.product(radii, concentrations)) do (R, Cs)
    L = U*T
    g(r) = snr(U, Cdiff, ∇Cdiff, r, R, C0, Cs, T, Dc, a) / ξ(U, r, T) - 1
    h = try
        find_zero(g, R, Order0())
    catch e
        R
    end
    h > R ? h : R
end

ic(C0) = IC.(30.0u"μm", radii.+a, S(C0).+a)

## figure
Rmin, Rmax = 0.5, 70.0
Cmin, Cmax = 2e-3, 1.0
Ic_str = rich("I", subscript("c"); font=:italic)
C_str = rich(rich("C", subscript("s"); font=:italic), " (μM)")
R_str = rich(rich("R"; font=:italic), " (μm)")
cmap = :viridis
clims = (0, 2.5)
clevels = range(clims...; step = 0.1)
fig = Figure(size=TwoColumns())
cb = Colorbar(fig[1,4],
    colormap = cmap,#cgrad(cmap, length(clevels)-1; categorical=true),
    colorrange = clims,
    ticks = 0:4,
    ticksvisible = false,
    tickformat = values -> [rich("10", superscript("$(Int(z))")) for z in values],
    label = Ic_str,
    #height = 0.85*h
)
C0 = [0, 10, 100]u"nM"
labs = rich.(["(A)", "(B)", "(C)"]; font=:bold)
ax = [
    Axis(fig[1,j], xlabel=R_str, ylabel=C_str, xscale=log10, yscale=log10,
        xticks=[1,3,9,27], yticks=[0.01,0.1,1],
        title=rich(labs[j], " C", subscript("0"), " = $(C0[j])")
    )
    for j in eachindex(C0)
]
for j in reverse(eachindex(C0))
    if j != 1
        ax[j].yticklabelsvisible = false
        ax[j].ylabelvisible = false
    end
    contourf!(ax[j],
        ustrip.(radii), ustrip.(concentrations), log10.(ic(C0[j]));
        levels=25, colormap=cmap, colorrange=clims
    )
    xlims!(ax[j], (Rmin, Rmax))
    ylims!(ax[j], (Cmin, Cmax))
end
fig
