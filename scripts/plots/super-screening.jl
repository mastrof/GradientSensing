using DrWatson
@quickactivate :GradientSensing
using DataFrames
using JLD2
using Roots
using GLMakie

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
ρ = 3.0u"μm"

S(ρ) = map(Iterators.product(radii, concentrations)) do (R, Cs)
    L = U*T
    f(r) = snr(U, Cexp, ∇Cexp, r*1u"μm", R, C0, Cs, ρ, T, Dc, a) - 1
    h = try
        #find_zero(f, R, Order0())
        find_zeros(f, ustrip(R), ustrip(40*R))[end] * 1u"μm"
    catch e
        R
    end
    h > R ? h : R
end

ic(ρ) = IC.(30.0u"μm", radii.+a, S(ρ).+a)
Ic = ic(ρ)

## Collect data for cuts along R and C
js = map(r -> findfirst(radii .≥ r), [1.0, 2.5, 5.0, 10.0, 25.0]u"μm")
Rx = radii[js]
I_R = Ic[js, :]

js = map(c -> findfirst(concentrations .≥ c), [0.015, 0.03, 0.05, 0.1]u"μM")
Cx = concentrations[js]
I_C = Ic[:, js]'


## figure
Rmin, Rmax = 0.5, 70.0
Cmin, Cmax = 2e-3, 1.0
Ic_str = rich("I", subscript("c"); font=:italic)
C_str = rich(rich("C", subscript("s"); font=:italic), " (μM)")
R_str = rich(rich("R"; font=:italic), " (μm)")
cmap = :viridis
clims = (0, 2.5)
clevels = range(clims...; step = 0.1)


## Initialize figure layout
fig = Figure(size = (1600,720))
pa = fig[1:2,1] = GridLayout()
pb = fig[1,2] = GridLayout()
pc = fig[2,2] = GridLayout()


## Panel B -- iso-R cuts
cmap = resample_cmap(:thermal, length(Rx))
labels = string.(round.(ustrip.(Rx), sigdigits=2))

ax_b = Axis(pb[1,1],
    xscale = log10,
    yscale = log10,
    xlabel = C_str,
    ylabel = Ic_str,
    xticks = [0.01, 0.1, 1],
    yticks = [1, 10, 100]
)
series!(ax_b,
    ustrip.(concentrations), I_R,
    color = cmap,
    linewidth = 9,
    labels = labels,
)
axislegend(ax_b; position=(0.02, 0),
    framevisible=true, framecolor=RGBAf(1.0,1.0,1.0,0.0),
    backgroundcolor=RGBAf(1.0,1.0,1.0,0.5)
)
text!(ax_b, 0.04, 0.8; text=R_str, space=:relative)
xlims!(ax_b, (Cmin, ustrip(concentrations[end])))
ylims!(ax_b, (1, 100))


## Panel C -- iso-C cuts
cmap = resample_cmap(:batlow, length(Cx))
labels = string.(round.(ustrip.(Cx), sigdigits=2))

ax_c = Axis(pc[1,1],
    xscale = log10,
    yscale = log10,
    xlabel = R_str,
    ylabel = Ic_str,
    xticks = [1, 3, 9, 27],
    yticks = [1, 10, 100]
)
series!(ax_c,
    ustrip.(radii), I_C,
    color = cmap,
    linewidth = 9,
    labels = labels
)
axislegend(ax_c; position=(1, 0.2))
text!(ax_c, 0.8, 0.72; text=C_str, space=:relative)
xlims!(ax_c, (Rmin, Rmax))
ylims!(ax_c, (1, 100))


## Panel A -- IC landscape
ax_a = Axis(pa[1,1],
    xlabel = R_str,
    ylabel = C_str,
    xscale = log10,
    yscale = log10,
    xticks = [1, 3, 9, 27],
    yticks = [0.01, 0.1, 1]
)

clims = (0, 2.5)
clevels = range(clims...; step = 0.1)
cmap = :viridis

cb = Colorbar(pa[1,2],
    colormap = cgrad(cmap, length(clevels)-1; categorical=true),
    colorrange = clims,
    ticks = 0:4,
    ticksvisible = false,
    tickformat = values -> [rich("10", superscript("$(Int(z))")) for z in values],
    label = Ic_str
)

contourf!(ax_a,
    ustrip.(radii), ustrip.(concentrations), log10.(Ic),
    colormap = cmap,
    levels = clevels,
    extendlow = :white
)
contour!(ax_a,
    ustrip.(radii), ustrip.(concentrations), log10.(Ic),
    levels = [clevels[findfirst(clevels .> 0)]],
    linewidth = 8,
    color = :white
)

scatter!(ax_a,
    ustrip.(Rx), repeat([1*Cmin], length(Rx)),
    marker = :utriangle,
    markersize = 45,
    strokewidth = 1,
    strokecolor = :white,
    color = resample_cmap(:thermal, length(Rx))
)

scatter!(ax_a,
    repeat([1*Rmin], length(Cx)), ustrip.(Cx),
    marker = :rtriangle,
    markersize = 45,
    strokewidth = 1,
    strokecolor = :white,
    color = resample_cmap(:batlow, length(Cx))
)

xlims!(ax_a, (Rmin, Rmax))
ylims!(ax_a, (Cmin, Cmax))

##
fig
