using DrWatson
@quickactivate :GradientSensing
using JLD2

function parameterspace_hein()
    # generate 100 log-uniform radii between 0.1 and 100 μm
    R = exp10.(range(-1, 2, length=100))u"μm"
    # use phytoplankton leakage rates to set bound values
    L₁ = leakage_rate(R[1], 0.01) |> u"pmol/s"
    L₂ = leakage_rate(R[end], 0.1) |> u"pmol/s"
    # generate 100 log-uniform values between L₁ and L₂
    L = exp10.(range(log10(ustrip(L₁)), log10(ustrip(L₂)), length=100))u"pmol/s"
    @save datadir("Hein", "paramspaceRL.jld2") R L
end

parameterspace_hein()