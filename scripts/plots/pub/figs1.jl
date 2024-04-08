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
Cexp(r, R, C0, Cs, ρ) = (C0 + Cs*R/r*exp(-(r-R)/ρ)) |> u"μM"
∇Cexp(r, R, C0, Cs, ρ) = (Cs*R*(ρ+r)*exp(-(r-R)/ρ) / (ρ*r^2)) |> u"μM/μm"
## sensing
snr(U, C, ∇C, r, R, C0, Cs, ρ, T, Dc, a) = signal(U,∇C,r,R,C0,Cs,ρ)/noise(U,C,r,R,C0,Cs,ρ,T,Dc,a)
signal(U, ∇C, r, R, C0, Cs, ρ) = U * ∇C(r,R,C0,Cs,ρ) |> u"μM/s"
noise(U, C, r, R, C0, Cs, ρ, T, Dc, a) = 6 * σ0(C,r,R,C0,Cs,ρ,T,Dc,a) * ξ(U,r,T) |> u"μM/s"
σ0(C,r,R,C0,Cs,ρ,T,Dc,a) = sqrt(3*C(r,R,C0,Cs,ρ) / (π*a*Dc*T^3*Unitful.Na))
#ξ(U, r, T) = sqrt(1 + 3/20*(U*T/r)^2)
ξ(U, r, T) = 1 / (1 - exp(-(r / (U*T))^1.5))

## parameter space
f = jldopen(datadir("HeinMod", "RC.jld2"), "r")
radii, concentrations = f["R"], f["Cₛ"]
close(f)

U = 50u"μm/s"
C0 = 1u"nM"
T = 100u"ms"
Dc = 500u"μm^2/s"
a = 0.5u"μm"

#=
SNR = let RR=[1.0, 2.5, 5.0, 20.0]u"μm", CC=[5.0, 50.0, 250.0, 1000.0]u"nM"
    map(Iterators.product(RR, CC)) do (R, Cs)
        r = range(R, 10R; length=1000)
        (r, @. snr(U, Cexp, ∇Cexp, r, R, C0, Cs, ρ, T, Dc, a))
    end
end
=#

S(ρ) = map(Iterators.product(radii, concentrations)) do (R, Cs)
    L = U*T
    f(r) = snr(U, Cexp, ∇Cexp, r, R, C0, Cs, ρ, T, Dc, a) - 1
    h = try
        find_zero(f, R)
    catch e
        R
    end
    h > R ? h : R
end

ic(ρ) = IC.(30.0u"μm", radii.+a, S(ρ).+a)

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
ax = [
    Axis(fig[1,j], xlabel=R_str, ylabel=C_str, xscale=log10, yscale=log10,
        xticks=[1,3,9,27], yticks=[0.01,0.1,1]
    )
    for j in 1:3
]
ρ = [30u"μm", 120u"μm", 400u"μm"]
labs = rich.(["(A)", "(B)", "(C)"]; font=:bold)
for j in 3:-1:1
    ax[j].title = rich(labs[j], " ", rich("ρ"; font=:italic), " = $(ustrip(ρ[j])) μm")
    if j != 1
        ax[j].yticklabelsvisible = false
        ax[j].ylabelvisible = false
    end
    contourf!(ax[j],
        ustrip.(radii), ustrip.(concentrations), log10.(ic(ρ[j]));
        levels=25, colormap=cmap, colorrange=clims
    )
    xlims!(ax[j], (Rmin, Rmax))
    ylims!(ax[j], (Cmin, Cmax))
end
fig
