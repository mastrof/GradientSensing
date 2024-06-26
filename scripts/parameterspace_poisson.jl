using DrWatson
@quickactivate :GradientSensing
using JLD2

"""
Defines the parameter space for the evaluation of sensing properties via direct
sampling of Poisson count events.

Returns a range of values for source radii (`R`) and source surface concentration (`Cₛ`),
spanning multiple orders of magnitude with a coarse spacing.
The values of `C` are chosen by first evaluating the leakage rates that would be
associated to a very small cell with low PER, and a large cell with medium-high PER,
and then evaluating the corresponding surface concentrations assuming
a background concentration of 1 nM (as defined in `src/global_constants.jl`);
the other values of `C` are selected by sampling uniformly in log scale between
the two extreme values thus obtained.

The values of `R` are saved in μm, those of `Cₛ` in μM.
"""
#== source radii ==#
# generate 20 log-uniform radii between 0.1 and 100 μm
R = exp10.(range(-1, 2, length=20))u"μm"
# add a few points in the low end
R = vcat(exp10.(range(-2,-1,length=5))u"μm", R)
# remove eventual duplicates and sort
unique!(R)
sort!(R)
# remove largest radius 100μm
R = R[1:end-1]

#== source concentration ==#
# fix some values for C₀ and Dc
C₀ = 1u"nM"
Dc = 608u"μm^2/s"
# get a small and a large value of phytoplankton leakage rates
# to set the bounds of the parameter space
R₁, R₂ = 0.1u"μm", 25u"μm"
Lmin = leakage_rate(R₁, 0.01) |> u"pmol/s"
Lmax = leakage_rate(R₂, 0.1) |> u"pmol/s"
# get equivalent Cmin and Cmax values
Cmin = C(R₁, Lmin, C₀, Dc)
Cmax = C(R₂, Lmax, C₀, Dc)
# generate 10 log-uniform values between Cmin and Cmax
Cₛ = exp10.(range(log10(ustrip(Cmin)), log10(ustrip(Cmax)), length=10))u"μM"
Cₛ = vcat(Cₛ, midpoints(Cₛ))
# remove eventual duplicates and sort
unique!(Cₛ)
sort!(Cₛ)

#== output ==#
@save datadir("Poisson", "RC.jld2") R Cₛ
