module GradientSensing

using Reexport
@reexport using Distributions, LinearAlgebra, Unitful, Statistics, StatsBase, Random
using HypothesisTests, Roots

Unitful.preferunits(u"s,μm,pmol"...)

include("global_constants.jl")
include("plankton_leakage.jl")
include("diffusive_field.jl")
include("sensing.jl")
include("phycosphere_Hein.jl")
include("chemotactic_index.jl")

end