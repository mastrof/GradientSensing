using DrWatson
@quickactivate :GradientSensing
using DataFrames
using JLD2
using Roots
using GLMakie

## concentration field
Cdiff(r, R, C0, Cs) = (C0 + Cs*R / r) |> u"μM"
∇Cdiff(r, R, C0, Cs) = (Cs*R / r^2) |> u"μM/μm"
## sensing
snr(U, C, ∇C, r, R, C0, Cs, T, Dc, a) = signal(U,∇C,r,R,C0,Cs)/noise(U,C,r,R,C0,Cs,T,Dc,a)
signal(U, ∇C, r, R, C0, Cs) = U * ∇C(r,R,C0,Cs) |> u"μM/s"
noise(U, C, r, R, C0, Cs, T, Dc, a) = 6 * σ0(C,r,R,C0,Cs,T,Dc,a) |> u"μM/s"
σ0(C,r,R,C0,Cs,T,Dc,a) = sqrt(3*C(r,R,C0,Cs) / (π*a*Dc*T^3*Unitful.Na))
ξ(U, r, T) = 1 / (1 - exp(-(r/(U*T))^(3/2)))

U = 50u"μm/s"
T = 100u"ms"
Dc = 500u"μm^2/s"
a = 0.5u"μm"
λ = 30.0u"μm"
C0 = 10u"nM"

Rs = [1.5, 8.0]u"μm"
Cs = [0.02, 0.1]u"μM"

## figure
fig = Figure(size=(1600,600), colormap=[:red, :red2, :cyan, :cyan2])
cmap = cgrad(:Paired_4, 4; categorical=true)
# Ic vs Dc
ax3 = Axis(fig[1,3]; yscale=log10, yticks=[1, 10, 100],
    xlabel=rich(rich("D", subscript("C"); font=:italic), " (μm", superscript("2"), "/s)"),
    yticklabelsvisible=false,
)
Ds = range(200, 1000; length=100)u"μm^2/s"
for (i, (C,R)) in enumerate(Iterators.product(Cs,Rs))
    ic = map(Ds) do Dc
        L = U*T
        g(r) = snr(U, Cdiff, ∇Cdiff, r, R, C0, C, T, Dc, a) / ξ(U, r, T) - 1
        h = try
            find_zero(g, R, Order0())
        catch e
            R
        end
        S = h > R ? h : R
        IC(λ, R+a, S+a)
    end
    lines!(ax3, collect(ustrip.(Ds)), ic, linewidth=10, color=cmap[i],
        label=rich(
            rich("R"; font=:italic), " = $(R)\n\n",
            rich("C", subscript("S"); font=:italic), " = $(C)",
        )
    )
end
# Ic vs U
ax2 = Axis(fig[1,2]; yscale=log10, yticks=[1, 10, 100],
    xlabel=rich(rich("U"; font=:italic), " (μm/s)"),
    yticklabelsvisible=false,
)
Us = range(15, 100; length=100)u"μm/s"
for (i, (C,R)) in enumerate(Iterators.product(Cs,Rs))
    ic = map(Us) do U
        L = U*T
        g(r) = snr(U, Cdiff, ∇Cdiff, r, R, C0, C, T, Dc, a) / ξ(U, r, T) - 1
        h = try
            find_zero(g, R, Order0())
        catch e
            R
        end
        S = h > R ? h : R
        IC(λ, R+a, S+a)
    end
    lines!(ax2, collect(ustrip.(Us)), ic, linewidth=10, color=cmap[i])
end
# Ic vs T
ax1 = Axis(fig[1,1]; yscale=log10, yticks=[1, 10, 100],
    xlabel=rich(rich("T"; font=:italic), " (ms)"),
    ylabel=rich("I", subscript("C"); font=:italic)
)
Ts = range(50, 500; length=100)u"ms"
for (i, (C,R)) in enumerate(Iterators.product(Cs,Rs))
    ic = map(Ts) do T
        L = U*T
        g(r) = snr(U, Cdiff, ∇Cdiff, r, R, C0, C, T, Dc, a) / ξ(U, r, T) - 1
        h = try
            find_zero(g, R, Order0())
        catch e
            R
        end
        S = h > R ? h : R
        IC(λ, R+a, S+a)
    end
    lines!(ax1, collect(ustrip.(Ts)), ic, linewidth=10, color=cmap[i])
end

Legend(fig[1,4], ax3; framevisible=false, rowgap=30)
linkyaxes!(ax1, ax2, ax3)
for (i,l) in enumerate(["A", "B", "C"])
    ax = i == 1 ? ax1 : (i == 2 ? ax2 : ax3)
    text!(ax, 0.05, 0.9; text=l, space=:relative, font=:bold)
end
fig
