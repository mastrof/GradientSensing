using DrWatson
@quickactivate :GradientSensing
using JLD2

"""
Evaluates phycosphere radii according to the Hein definition (SNR=1)
for all values of `R` and `Cₛ` defined in `paramspaceRC.jld2`

Returns a jld2 file with a phycosphere radius (`S`) for each `(R,Cₛ)` pair.
The data produced herein is then used to evaluate chemotactic indices (`scripts/ic_hein.jl`)
"""
function phycosphere_hein()
    f = jldopen(datadir("Hein", "paramspaceRC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]
    
    # set T, U, Π to some default values
    T = 100u"ms"
    U = 46.5u"μm/s"
    Π = 6u"1"
    params = @strdict C₀ T U Π
    # need to have a ustrip copy for the savename
    params_ustrip = Dict(keys(params) .=> ustrip.(values(params)))
    
    produce_or_load(
        datadir("Hein"), params_ustrip;
        prefix="phycosphere", suffix="jld2",
        tag = false
    ) do params_ustrip
        # evaluate phycosphere radius for each (R,L) pair
        S = map(p -> HeinRadius(p..., C₀, T, U, Π), Iterators.product(R,Cₛ))
        @strdict S
    end
end

phycosphere_hein()