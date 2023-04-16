module GradientSensing

using Reexport
@reexport using Distributions, LinearAlgebra, Unitful, Statistics, StatsBase, Random
using HypothesisTests, Roots
using Plots; Plots.pyplot()

Unitful.preferunits(u"s,μm,pmol"...)

Plots.default(
    grid=false, frame=:border, thickness_scaling=1.5,
    bgcolorlegend=:transparent, fgcolorlegend=:transparent,
    palette=:Dark2
)

include("global_constants.jl")
include("plankton_leakage.jl")
include("diffusive_field.jl")
include("sensing.jl")
include("phycosphere_Hein.jl")
include("chemotactic_index.jl")

end