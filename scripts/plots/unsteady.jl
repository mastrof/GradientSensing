using DrWatson
@quickactivate :GradientSensing
using DataFrames
using JLD2
using Roots
using SpecialFunctions
using GLMakie

## concentration field
function Cexp(r, R, C0, Cs, D, t)
    Cst = C0 + Cs*R/r
    unst = erfc((r-R) / (2*sqrt(D*t)))
    (Cst * unst) |> u"μM"
end
function ∇Cexp(r, R, C0, Cs, D, t)
    Cst = C0 + Cs*R/r
    unst = erfc((r-R) / (2*sqrt(D*t)))
    ∇Cst = Cs*R/r^2
    unst_dt = 1/sqrt(π*D*t) * exp(-(r-R)^2 / (4*D*t))
    (∇Cst*unst + Cst*unst_dt) |> u"μM/μm"
end
## sensing
snr(U, C, ∇C, r, R, C0, Cs, T, Dc, a, t) = signal(U,∇C,r,R,C0,Cs,Dc,t)/noise(U,C,r,R,C0,Cs,T,Dc,a,t)
signal(U, ∇C, r, R, C0, Cs, Dc, t) = U * ∇C(r,R,C0,Cs,Dc,t) |> u"μM/s"
noise(U, C, r, R, C0, Cs, T, Dc, a, t) = 6 * σ0(C,r,R,C0,Cs,T,Dc,a,t) * ξ(U,r,T) |> u"μM/s"
σ0(C,r,R,C0,Cs,T,Dc,a,t) = sqrt(3*C(r,R,C0,Cs,Dc,t) / (π*a*Dc*T^3*Unitful.Na))
#ξ(U, r, T) = sqrt(1 + 3/20*(U*T/r)^2)
ξ(U, r, T) = 1 / (1 - exp(-(r / (U*T))^1.5))

## parameter space
f = jldopen(datadir("HeinMod", "RC.jld2"), "r")
radii, concentrations = f["R"][1:2:end], f["Cₛ"][1:2:end]
close(f)

U = 50u"μm/s"
C0 = 1u"nM"
T = 100u"ms"
Dc = 500u"μm^2/s"
a = 0.5u"μm"

# converge to steady state
let
    r = radii
    t = range(1, 600; length=1000)u"s"
    τ = @. erfc(r/(2*sqrt(Dc*t')))
    fig = Figure(size=(800, 600))
    xlabel = rich(rich("r-R"; font=:italic), " (μm)")
    ylabel = rich(rich("t"; font=:italic), " (s)")
    cblabel = rich("C(r,t) / C(r,t→∞)"; font=:italic)
    ax = Axis(fig[1,1];
        xlabel, ylabel,
        yscale=log10, yticks=[1, 10, 100],
    )
    xlims!(ustrip.(extrema(r))...)
    ylims!(ustrip.(extrema(t))...)
    plt = contourf!(ax, ustrip(r), ustrip(t), ustrip(τ))
    contour!(ax, ustrip(r), ustrip(t), ustrip(τ);
        levels=[0.9], color=:black, linewidth=20
    )
    lines!(ax, ustrip(r), r -> 125(r^2/(4*ustrip(Dc)));
        color=:white, linewidth=10
    )
    Colorbar(fig[1,2], plt; label=cblabel, ticks=0:0.2:0.8)
    fig
end


# unsteady ic
S(n) = map(Iterators.product(radii, concentrations)) do (R, Cs)
    #τ = (R^2 / 2Dc) |> u"s"
    τ = 1u"s"
    #f(r) = snr(U, Cexp, ∇Cexp, r*1u"μm", R, C0, Cs, T, Dc, a, n*τ) - 1
    f(r) = snr(U, Cexp, ∇Cexp, r, R, C0, Cs, T, Dc, a, n*τ) - 1
    h = try
        find_zero(f, 2R, Order0())
        #z = find_zeros(f, ustrip(R), ustrip(50*R))
        #z[end] * 1u"μm"
    catch e
        #println(e)
        R
    end
    h > R ? h : R
end

ic(n) = IC.(30.0u"μm", radii.+a, S(n).+a)

## figure
Rmin, Rmax = (0.5, 70.0)
Cmin, Cmax = (2e-3, 1.0)
Ic_str = rich("I", subscript("c"); font=:italic)
C_str = rich(rich("C", subscript("s"); font=:italic), " (μM)")
R_str = rich(rich("R"; font=:italic), " (μm)")
cmap = :viridis
#clims = (0, 2.5)
clims = (0, 100)
clevels = range(clims...; step = 0.1)
fig = Figure(size=(1600,600))
cb = Colorbar(fig[1,4],
    colormap = cmap,#cgrad(cmap, length(clevels)-1; categorical=true),
    colorrange = clims,
    #ticks = 0:4,
    ticksvisible = false,
    #tickformat = values -> [rich("10", superscript("$(Int(z))")) for z in values],
    label = Ic_str,
)
ns = [1, 5, Inf]
ax = [
    Axis(fig[1,n];
        xlabel=R_str, ylabel=C_str, xscale=log10, yscale=log10,
        xticks=[1,3,9,27], yticks=[0.01,0.1,1],
        title=rich(rich("t"; font=:bold_italic), " = $(ns[n]) s"),
    )
    for n in eachindex(ns)
]
for i in eachindex(ax)
    if i != 1
        ax[i].ylabel = ""
    end
    n = ns[i]
    contourf!(ax[i],
        ustrip.(radii), ustrip.(concentrations), ustrip(S(n) .- radii);
        levels=range(0,100;length=10)
        #levels=25, colormap=cmap, colorrange=clims
        #colormap=cmap, levels=clevels
    )
    xlims!(ax[i], (Rmin, Rmax))
    ylims!(ax[i], (Cmin, Cmax))
end
fig
