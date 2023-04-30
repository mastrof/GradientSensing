module GradientSensing

using Reexport
@reexport using Distributions, LinearAlgebra, Unitful, Statistics, StatsBase, Random
using HypothesisTests, Roots

# import unitful dimension types for dispatch
using Unitful: 𝐍, 𝐋, 𝐓
# set preferred units
Unitful.preferunits(u"s,μm,pmol"...)

include("global_constants.jl")
include("plankton_leakage.jl")
include("diffusive_field.jl")
include("sensing.jl")
include("phycosphere_Hein.jl")

include("chemotactic_index.jl")

include("poisson_events.jl")

end