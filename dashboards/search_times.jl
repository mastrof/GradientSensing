using DrWatson
@quickactivate "GradientSensing"
using GradientSensing
using JLD2
using DataFrames
using Roots
using Interpolations
using StatsBase
using GLMakie
set_theme!(Theme(
    fontsize = 32,
    palette = (
        color = cgrad(:Dark2_8; categorical=true),
        marker = [:circle,:rect,:diamond,:utriangle,:dtriangle,:rtriangle,:ltriangle,:hexagon,:pentagon,:cross,:xcross,:star4,:star5,:star6,:star8],
        linestyle = :solid,
    ),
    Axis = (
        xgridvisible = false, ygridvisible = false,
        xticksize = -10, yticksize = -10,
        xminorticksvisible = false, yminorticksvisible = false,
        #xticksmirrored = true, yticksmirrored = true,
        titlefont = :bold, titlesize = 32,
    ),
    Legend = (
        framevisible = false,
        titlefont = :regular,
        titlesize = 28, labelsize = 28,
    ),
    Colorbar = (
        ticksvisible=false,
    ),
    Label = (
        font = :bold, fontsize = 24, halign = :center,
    ),
    Lines = (
        linewidth = 4, cycle = Cycle([:color]),
    ),
    Scatter = (
        markersize = 16, cycle = Cycle([:color, :marker], covary=true),
    ),
    ScatterLines = (
        linewidth = 4, markersize = 16,
        cycle = Cycle([:color, :marker], covary=true),
    ),
))


# concentration field
Cfield(r, R, C0, Cs) = (C0 + Cs*R / r) |> u"μM"
Cgrad(r, R, C0, Cs) = (Cs*R / r^2) |> u"μM/μm"
# sensing
snr(U,C,∇C,r,R,C0,Cs,T,Dc,a) = signal(U,∇C,r,R,C0,Cs)/noise(U,C,r,R,C0,Cs,T,Dc,a)
signal(U,∇C,r,R,C0,Cs) = U * ∇C(r,R,C0,Cs) |> u"μM/s"
noise(U,C,r,R,C0,Cs,T,Dc,a) = 6*σ0(C,r,R,C0,Cs,T,Dc,a) |> u"μM/s"
σ0(C,r,R,C0,Cs,T,Dc,a) = sqrt(3*C(r,R,C0,Cs) / (π*a*Dc*T^3*Unitful.Na))
ξ(U,r,T) = 1 / (1 - exp(-(r/(U*T))^(3/2)))
# phytoplankton
function community_structure(α, Ntot, nclasses, Rmin, Rmax)
    V_min, V_max = @. 4π / 3 * (Rmin^3, Rmax^3) # μm^3
    V_edges = exp2.(range(log2(V_min), log2(V_max); length=nclasses+1))
    V = midpoints(V_edges)
    r_edges = @. (3 * V_edges / (4π))^(1/3)
    r = @. (3 * V / (4π))^(1/3) # μm
    N = @. V^α
    N = N ./ sum(N) .* Ntot # cells / mL
    return r, r_edges, N
end
# search
function encounter_time(N, R, λ, U, T, Dc, C0, a, PER; Π=6)
    Te = time_to_encounter(ustrip(R), ustrip(λ), ustrip(U), N) / 3600 # hours
    L = leakage_rate(R, PER)
    Cs = C(R, L, C0, Dc) - C0
    S = HeinModRadius(R, Cs, C0, T, Dc, U, Π) |> u"μm"
    Ic = IC(λ, R+a, S+a)
    Te / Ic
end
function time_to_encounter(R, λ, U, N)
    Γ = encounter_kernel(R, λ, U) # μm^3/s
    Γ *= 1e-12 # mL/s
    Te = 1 / (Γ*N) # s
    return Te
end
function encounter_kernel(R, λ, U)
    kn = 4λ / 3R
    β = kn / (1 + kn)
    return β * π * U * R^2
end

# global parameters
nclasses = 15
Rmin, Rmax = 0.5, 70.0
PER = [0.02, 0.1, 0.4]

# figure setup
fig = Figure(size=(1600,900))
# 12 rows, 12 columns
pa = fig[3:10, 2:5] = GridLayout() # community structure
# use 5:8, 1:2 for phytoplankton sliders
pb = fig[3:10, 6:12] = GridLayout() # search times

# set parameters with sliders + labels
## chemotaxis - global ; put on top
sl_Dc = Slider(fig[12,3:4]; range=range(100,1000,step=100), startvalue=500, width=200)
Dc = @lift($(sl_Dc.value) * 1u"μm^2/s")
lab_Dc = Label(fig[11,3:4]; text=@lift("Dc = $($Dc)"), padding=(0,0,0,25))
## chemotaxis - 1
c1 = cgrad(:Dark2)[1]
sl_U1 = Slider(fig[2,6]; range=range(15,100,step=5), startvalue=50, width=200)
sl_T1 = Slider(fig[2,8]; range=range(50,500,step=50), startvalue=100, width=200)
sl_a1 = Slider(fig[2,10]; range=range(0.25,2.5,step=0.25), startvalue=0.5, width=200)
sl_λ1 = Slider(fig[2,12]; range=range(15,100,step=5), startvalue=30, width=200)
U1 = @lift($(sl_U1.value) * 1u"μm/s")
T1 = @lift($(sl_T1.value) * 1u"ms")
a1 = @lift($(sl_a1.value) * 1u"μm")
λ1 = @lift($(sl_λ1.value) * 1u"μm")
lab_U1 = Label(fig[1,6]; text=@lift("U = $($U1)"), padding=(0,0,0,25), color=c1)
lab_T1 = Label(fig[1,8]; text=@lift("T = $($T1)"), padding=(0,0,0,25), color=c1)
lab_a1 = Label(fig[1,10]; text=@lift("a = $($a1)"), padding=(0,0,0,25), color=c1)
lab_λ1 = Label(fig[1,12]; text=@lift("λ = $($λ1)"), padding=(0,0,0,25), color=c1)
## chemotaxis - 2
c2 = cgrad(:Dark2)[3]
sl_U2 = Slider(fig[12,6]; range=range(15,100,step=5), startvalue=50, width=200)
sl_T2 = Slider(fig[12,8]; range=range(50,500,step=50), startvalue=100, width=200)
sl_a2 = Slider(fig[12,10]; range=range(0.25,2.5,step=0.25), startvalue=0.5, width=200)
sl_λ2 = Slider(fig[12,12]; range=range(15,100,step=5), startvalue=30, width=200)
U2 = @lift($(sl_U2.value) * 1u"μm/s")
T2 = @lift($(sl_T2.value) * 1u"ms")
a2 = @lift($(sl_a2.value) * 1u"μm")
λ2 = @lift($(sl_λ2.value) * 1u"μm")
lab_U2 = Label(fig[11,6]; text=@lift("U = $($U2)"), padding=(0,0,0,25), color=c2)
lab_T2 = Label(fig[11,8]; text=@lift("T = $($T2)"), padding=(0,0,0,25), color=c2)
lab_a2 = Label(fig[11,10]; text=@lift("a = $($a2)"), padding=(0,0,0,25), color=c2)
lab_λ2 = Label(fig[11,12]; text=@lift("λ = $($λ2)"), padding=(0,0,0,25), color=c2)
## phytoplankton - oligotrophic
cp1 = cgrad(:tab10)[1]
sl_N1 = Slider(fig[2,1]; range=range(2e4,1e6,step=2e4), startvalue=6e4, width=200)
sl_α1 = Slider(fig[4,1]; range=range(0.65,1.2,step=0.05), startvalue=1.0, width=200)
sl_C01 = Slider(fig[6,1]; range=range(0,500,step=10), startvalue=10, width=200)
N1 = @lift($(sl_N1.value) * 1)
α1 = @lift($(sl_α1.value) * 1)
C01 = @lift($(sl_C01.value) * 1u"nM")
lab_N1 = Label(fig[1,1]; text=@lift("N = $($N1) cells/mL"), padding=(0,0,0,25), color=cp1)
lab_α1 = Label(fig[3,1]; text=@lift("α = $($α1)"), padding=(0,0,0,25), color=cp1)
lab_C01 = Label(fig[5,1]; text=@lift("C₀ = $($C01)"), padding=(0,0,0,25), color=cp1)
## phytoplankton - productive
cp2 = cgrad(:tab10)[2]
sl_N2 = Slider(fig[8,1]; range=range(2e4,1e6,step=2e4), startvalue=2e5, width=200)
sl_α2 = Slider(fig[10,1]; range=range(0.65,1.2,step=0.05), startvalue=0.75, width=200)
sl_C02 = Slider(fig[12,1]; range=range(0,500,step=10), startvalue=100, width=200)
N2 = @lift($(sl_N2.value) * 1)
α2 = @lift($(sl_α2.value) * 1)
C02 = @lift($(sl_C02.value) * 1u"nM")
lab_N2 = Label(fig[7,1]; text=@lift("N = $($N2) cells/mL"), padding=(0,0,0,25), color=cp2)
lab_α2 = Label(fig[9,1]; text=@lift("α = $($α2)"), padding=(0,0,0,25), color=cp2)
lab_C02 = Label(fig[11,1]; text=@lift("C₀ = $($C02)"), padding=(0,0,0,25), color=cp2)

# evaluate community structure
comm1 = @lift(community_structure(-$α1, $N1, nclasses, Rmin, Rmax))
r1 = @lift($(comm1)[1] .* 1u"μm")
redg1 = @lift($(comm1)[2] .* 1u"μm")
abundance1 = @lift($(comm1)[3])
comm2 = @lift(community_structure(-$α2, $N2, nclasses, Rmin, Rmax))
r2 = @lift($(comm2)[1] .* 1u"μm")
redg2 = @lift($(comm2)[2] .* 1u"μm")
abundance2 = @lift($(comm2)[3])

# panel A - community structure
ax_A = Axis(pa[1,1];
    xlabel=rich("Phytoplankton radius ", rich("R"; font=:italic), " (μm)"),
    ylabel=rich("Abundance ", rich("N"; font=:italic), " (cells/mL)"),
    xscale=log10, yscale=log10,
    xticks=[1, 3, 9, 27], yticks=exp10.(1:2:6),
    xtickformat=(xs -> string.(Int.(xs))),
    ytickformat=(ys -> [rich("10", superscript("$(Int(log10(y)))")) for y in ys]),
    title="Community structure"
)
ylims!(ax_A, (1, 5e5))
xlims!(ax_A, (Rmin, Rmax))
scatterlines!(ax_A, ustrip.(r1[]), abundance1;
    linewidth=8, marker=:utriangle, markersize=32, color=cp1,
    label="Oligotrophic"
)
scatterlines!(ax_A, ustrip.(r2[]), abundance2;
    linewidth=8, marker=:dtriangle, markersize=32, color=cp2,
    label="Productive"
)
axislegend(ax_A; position=:rt, patchsize=(60, 20))

# evaluate search times
Te_rnd_oli1 = @lift(
    @. encounter_time($abundance1, $r1, $λ1, $U1, $T1, $Dc, 0u"nM", $a1, 0)
)
Te_rnd_oli2 = @lift(
    @. encounter_time($abundance1, $r1, $λ2, $U2, $T2, $Dc, 0u"nM", $a2, 0)
)
Te_rnd_pro1 = @lift(
    @. encounter_time($abundance2, $r2, $λ1, $U1, $T1, $Dc, 0u"nM", $a1, 0)
)
Te_rnd_pro2 = @lift(
    @. encounter_time($abundance2, $r2, $λ2, $U2, $T2, $Dc, 0u"nM", $a2, 0)
)
Te_oli1 = [
    @lift(@. encounter_time($abundance1, $r1, $λ1, $U1, $T1, $Dc, $C01, $a1, $(PER)[i]))
    for i in eachindex(PER)
]
Te_oli2 = [
    @lift(@. encounter_time($abundance1, $r1, $λ2, $U2, $T2, $Dc, $C01, $a2, $(PER)[i]))
    for i in eachindex(PER)
]
Te_pro1 = [
    @lift(@. encounter_time($abundance2, $r2, $λ1, $U1, $T1, $Dc, $C02, $a1, $(PER)[i]))
    for i in eachindex(PER)
]
Te_pro2 = [
    @lift(@. encounter_time($abundance2, $r2, $λ2, $U2, $T2, $Dc, $C02, $a2, $(PER)[i]))
    for i in eachindex(PER)
]

# plot search times in oligotrophic waters (environment 1)
ax_oli = Axis(pb[1,1];
    xlabel=rich(rich("R"; font=:italic), " (μm)"),
    ylabel=rich(rich("T", subscript("e"); font=:italic), " (hours)"),
    yscale=log10,
    xscale=log10,
    yticks=exp10.(0:3),
    xticks=[1, 3, 9, 27],
    title="Oligotrophic waters",
    titlecolor=cp1
)
xlims!(ax_oli, (0.45, 75))
ylims!(ax_oli, (1, 2e3))
lines!(ax_oli, ustrip.(r1[]), Te_rnd_oli1;
    linestyle=:dash, linewidth=3, color=c1
)
lines!(ax_oli, ustrip.(r1[]), Te_rnd_oli2;
    linestyle=:dash, linewidth=3, color=c2
)
for (i,Te) in enumerate((Te_oli1, Te_oli2))
    color = i == 1 ? c1 : c2
    band!(ax_oli, ustrip.(r1[]), Te[1], Te[3]; color=color, alpha=0.3)
    scatterlines!(ax_oli, ustrip.(r1[]), Te[2];
        color=color, linewidth=8, markersize=24,
        marker=(i == 1 ? :circle : :rect)
    )
end

# plot search times in productive waters (environment 2)
ax_pro = Axis(pb[1,2];
    xlabel=rich(rich("R"; font=:italic), " (μm)"),
    #ylabel=rich(rich("T", subscript("e"); font=:italic), " (hours)"),
    yscale=log10,
    xscale=log10,
    yticks=exp10.(0:3), yticklabelsvisible=false,
    xticks=[1, 3, 9, 27],
    title="Productive waters",
    titlecolor=cp2
)
xlims!(ax_pro, (0.45, 75))
ylims!(ax_pro, (1, 2e3))
lines!(ax_pro, ustrip.(r2[]), Te_rnd_pro1;
    linestyle=:dash, linewidth=3, color=c1
)
lines!(ax_pro, ustrip.(r2[]), Te_rnd_pro2;
    linestyle=:dash, linewidth=3, color=c2
)
for (i,Te) in enumerate((Te_pro1, Te_pro2))
    color = i == 1 ? c1 : c2
    band!(ax_pro, ustrip.(r2[]), Te[1], Te[3]; color=color, alpha=0.3)
    scatterlines!(ax_pro, ustrip.(r2[]), Te[2];
        color=color, linewidth=8, markersize=24,
        marker=(i == 1 ? :circle : :rect)
    )
end
