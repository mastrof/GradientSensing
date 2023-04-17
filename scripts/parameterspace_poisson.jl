using DrWatson
@quickactivate :GradientSensing
using JLD2

function parameterspace_poisson()
    # generate 20 log-uniform radii between 0.5 and 100 μm
    R = exp10.(range(log10(0.5), 2, length=20))u"μm"
    # use phytoplankton leakage rates to set bound values
    L₁ = leakage_rate(R[1], 0.01) |> u"pmol/s"
    L₂ = leakage_rate(R[end], 0.01) |> u"pmol/s"
    # generate 10 log-uniform values between L₁ and L₂
    L = exp10.(range(log10(ustrip(L₁)), log10(ustrip(L₂)), length=10))u"pmol/s"
    @save datadir("PoissonSampling", "paramspaceRL.jld2") R L
end

parameterspace_poisson()