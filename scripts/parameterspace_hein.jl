using DrWatson
@quickactivate :GradientSensing
using JLD2

"""
Defines the parameter space for the evaluation of sensing properties using
the Hein definition of phycosphere (SNR=1).

Returns a range of values for source radii (`R`) and source surface concentration (`Cₛ`),
spanning multiple orders of magnitude.
The values of `C` are chosen by first evaluating the leakage rates that would be
associated to a very small cell with low PER, and a large cell with medium-high PER,
and then evaluating the corresponding surface concentrations assuming
a background concentration of 1 nM and a molecular diffusivity of 608 μm^2/s;
the other values of `C` are selected by sampling uniformly in log scale between
the two extreme values thus obtained.

The values of `R` are saved in μm, those of `Cₛ` in μM.
"""




#== source radii ==#
# generate 100 log-uniform radii between 0.1 and 100 μm
R = exp10.(range(-1, 2, length=100))u"μm"

#== source concentration ==#
# get a small and a large value of phytoplankton leakage rates
# to set the bounds of the parameter space
R₁, R₂ = 0.1u"μm", 100u"μm"
PER₁, PER₂ = 0.01, 0.1
Lmin = leakage_rate(R₁, PER₁) |> u"pmol/s"
Lmax = leakage_rate(R₂, PER₂) |> u"pmol/s"
# get equivalent Cmin and Cmax values using standard values for C₀ and Dc
C₀ = 1u"nM"
Dc = 608u"μm^2/s"
Cmin = C(R₁, Lmin, C₀, Dc)
Cmax = C(R₂, Lmax, C₀, Dc)
# generate 100 log-uniform values in this range
Cₛ = exp10.(range(log10(ustrip(Cmin)), log10(ustrip(Cmax)), length=100))u"μM"

@save datadir("HeinMod", "RC.jld2") R Cₛ
