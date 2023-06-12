using DrWatson
@quickactivate :GradientSensing
using JLD2

"""
Evaluates phycosphere radii according to the Hein definition (SNR=1),
modulated by a custom function to account for Nyquist frequency,
for all values of `R` and `Cₛ` defined in `RC.jld2`

Returns a jld2 file with a phycosphere radius (`S`) for each `(R,Cₛ)` pair.
"""
function phycosphere_hein()
    f = jldopen(datadir("HeinMod", "RC.jld2"))
    R, Cₛ = f["R"], f["Cₛ"]

    # set T, U, Π to some default values
    T = 100u"ms"
    U = 46.5u"μm/s"
    Π = 6u"1"
    params = @strdict C₀ T U Π
    # need to have a ustrip copy for the savename
    params_ustrip = Dict(keys(params) .=> ustrip.(values(params)))

    produce_or_load(
        datadir("HeinMod"), params_ustrip;
        prefix="phycosphere", suffix="jld2",
        tag = false, force = true
    ) do params_ustrip
        # evaluate phycosphere radius for each (R,L) pair
        S = map(p -> HeinModRadius(p..., C₀, T, U, Π), Iterators.product(R,Cₛ))
        @strdict S
    end
end

phycosphere_hein()
