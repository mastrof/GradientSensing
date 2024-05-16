using DrWatson
@quickactivate "GradientSensing"
begin
	using GradientSensing
	using JLD2
	using DataFrames
	using Roots
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
end

begin
	# concentration field
	Cfield(r, R, C0, Cs) = (C0 + Cs*R / r) |> u"μM"
	Cgrad(r, R, C0, Cs) = (Cs*R / r^2) |> u"μM/μm"
	# sensing
	snr(U,C,∇C,r,R,C0,Cs,T,Dc,a) = signal(U,∇C,r,R,C0,Cs)/noise(U,C,r,R,C0,Cs,T,Dc,a)
	signal(U,∇C,r,R,C0,Cs) = U * ∇C(r,R,C0,Cs) |> u"μM/s"
	noise(U,C,r,R,C0,Cs,T,Dc,a) = 6*σ0(C,r,R,C0,Cs,T,Dc,a) |> u"μM/s"
	σ0(C,r,R,C0,Cs,T,Dc,a) = sqrt(3*C(r,R,C0,Cs) / (π*a*Dc*T^3*Unitful.Na))
	ξ(U,r,T) = 1 / (1 - exp(-(r/(U*T))^(3/2)))
end

begin
	# define parameter space
	Rs = exp10.(range(-1, 2, length=40))u"μm"
	Cs = exp10.(range(-2, 0, length=40))u"μM"
end

begin
	# figure setup
	fig = Figure(size=(1600,900))
    # 12 rows needed for sliders + labels
	pa = fig[1:12, 1:3] = GridLayout()
	pb = fig[1:6, 4:5] = GridLayout()
	pc = fig[7:12, 4:5] = GridLayout()
	
	# global parameters
	Rmin, Rmax = 0.5, 70.0
	Cmin, Cmax = 1e-2, 1.0
	clims = (0, 2.5)
	clevels = range(clims...; step=0.1)
	Ic_str = rich("I", subscript("C"); font=:italic)
	C_str = rich("Excess concentration ", rich("C", subscript("S"); font=:italic), " (μM)")
	R_str = rich("Phytoplankton radius ", rich("R"; font=:italic), " (μm)")

	# set parameters with sliders
	sl_U = Slider(fig[2,0]; range=range(15,100,step=5), startvalue=50, width=200)
	sl_T = Slider(fig[4,0]; range=range(50,500,step=50), startvalue=100)
	sl_Dc = Slider(fig[6,0]; range=range(100,1000,step=100), startvalue=500)
	sl_a = Slider(fig[8,0]; range=range(0.25,2.5,step=0.25), startvalue=0.5)
	sl_C0 = Slider(fig[10,0]; range=range(0,500,step=10), startvalue=100)
	# ideally λ must be larger than UT...
	sl_λ = Slider(fig[12,0]; range=range(15,100,step=5), startvalue=30)
	# define variables and assign units
	U = @lift($(sl_U.value) * 1u"μm/s")
	T = @lift($(sl_T.value) * 1u"ms")
	Dc = @lift($(sl_Dc.value) * 1u"μm^2/s")
	a = @lift($(sl_a.value) * 1u"μm")
	C0 = @lift($(sl_C0.value) * 1u"nM")
	λ = @lift($(sl_λ.value) * 1u"μm")
	L = @lift($U * $T) # sensory lengthscale

    # slider labels
    lab_U = Label(fig[1,0]; text=@lift("U = $($U)"), padding=(0,0,0,25))
    lab_T = Label(fig[3,0]; text=@lift("T = $($T)"), padding=(0,0,0,25))
    lab_Dc = Label(fig[5,0]; text=@lift("Dc = $($Dc)"), padding=(0,0,0,25))
    lab_a = Label(fig[7,0]; text=@lift("a = $($a)"), padding=(0,0,0,25))
    lab_C0 = Label(fig[9,0]; text=@lift("C₀ = $($C0)"), padding=(0,0,0,25))
    lab_λ = Label(fig[11,0]; text=@lift("λ = $($λ)"), padding=(0,0,0,25))

	# compute IC
    g(r,U,R,C0,C,T,Dc,a) = snr(U,Cfield,Cgrad,r,R,C0,C,T,Dc,a) / ξ(U,r,T) - 1
    ic = Observable(ones(length(Rs), length(Cs)))
    ic_log10 = @lift(log10.($ic))
    RCSpace = Iterators.product(Rs, Cs)
    observer = @lift([$U,$T,$Dc,$a,$C0,$λ])
    on(observer) do _
        for (i,(R,C)) in enumerate(RCSpace)
            ic[][i] = begin
                h = try
                    find_zero(r -> g(r,U[],R,C0[],C,T[],Dc[],a[]), R, Order5())
                catch e
                    R
                end
                S = h > R ? h : R
                # evaluate IC
                IC(λ[],R+a[],S+a[]) # defined in GradientSensing
            end
            ic_log10[][i] = log10(ic[][i])
        end
        notify(ic)
        notify(ic_log10)
    end

    # plot cuts at selected values of R
    ax_B = Axis(pb[1,1];
        xlabel=rich(C_str.children[2:end]...), ylabel=Ic_str,
        xscale=log10, yscale=log10,
        xticks=[0.01, 0.1, 1],
        yticks=[1, 10, 100]
    )
    js = [11, 18, 25, 32]
    R_cuts = Rs[js]
    R_cut_labels = @. string(round(typeof(1.0u"μm"), R_cuts; sigdigits=1))
    ic_cuts_R = [@lift(view($(ic), j, :)) for j in js]
    cmap = cgrad(:thermal, length(R_cuts); categorical=true)
    for j in eachindex(R_cuts)
        lines!(ax_B, ustrip.(Cs), ic_cuts_R[j];
            linewidth=8, color=cmap[j], label=R_cut_labels[j]
        )
    end
    axislegend(ax_B; position=(0.02, 1))
    xlims!(ax_B, (Cmin, Cmax))
    ylims!(ax_B, (1, 100))

    # plot cuts at selected values of C
    ax_C = Axis(pc[1,1];
        xlabel=rich(R_str.children[2:end]...), ylabel=Ic_str,
        xscale=log10, yscale=log10,
        xticks=[1, 3, 9, 27],
        yticks=[1, 10, 100]
    )
    ks = [5, 11, 17, 23]
    C_cuts = Cs[ks]
    C_cut_labels = @. string(round(typeof(1.0u"μM"), C_cuts; sigdigits=1))
    ic_cuts_C = [@lift(view($(ic), :, k)) for k in ks]
    cmap = cgrad(:batlow, length(C_cuts); categorical=true)
    for k in eachindex(C_cuts)
        lines!(ax_C, ustrip.(Rs), ic_cuts_C[k];
            linewidth=8, color=cmap[k], label=C_cut_labels[k]
        )
    end
    axislegend(ax_C; position=(1, 1))
    xlims!(ax_C, (Rmin, Rmax))
    ylims!(ax_C, (1, 100))
	
	# plot IC landscape
	cmap = :viridis
	ax_A = Axis(pa[1,1];
		xlabel=R_str, ylabel=C_str,
		xscale=log10, yscale=log10,
		xticks=[1, 3, 9, 27],
		yticks=[0.01, 0.1, 1],
	)
	xlims!(ax_A, (Rmin, Rmax))
	ylims!(ax_A, (Cmin, Cmax))
	contourf!(ax_A, ustrip.(Rs), ustrip.(Cs), ic_log10;
		colormap=cmap, levels=clevels, extendlow=:white, extendhigh=:red,
	)
	cb = Colorbar(pa[1,2];
		colormap=cgrad(cmap, length(clevels)-1; categorical=true),
		colorrange=clims, highclip=:red,
		ticks=0:4, ticksvisible=false,
		tickformat = zs -> [rich("10", superscript("$(Int(z))")) for z in zs],
		label = Ic_str
	)

    # cut markers
    scatter!(ax_A, ustrip.(R_cuts), repeat([Cmin], length(R_cuts));
        marker=:utriangle, markersize=45,
        strokewidth=1, strokecolor=:white,
        color=cgrad(:thermal, length(R_cuts); categorical=true).colors.colors
    )
    scatter!(ax_A, repeat([Rmin], length(C_cuts)), ustrip.(C_cuts);
        marker=:rtriangle, markersize=45,
        strokewidth=1, strokecolor=:white,
        color=cgrad(:batlow, length(C_cuts); categorical=true).colors.colors
    )

	fig
end
