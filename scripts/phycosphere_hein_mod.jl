##
using DrWatson
@quickactivate :GradientSensing
using JLD2

##
"""
Evaluates phycosphere radii according to the Hein definition (SNR=1),
modulated by a custom function to account for Nyquist frequency,
for all values of `R` and `Cₛ` defined in `RC.jld2`

Returns a jld2 file with a phycosphere radius (`S`) for each `(R,Cₛ)` pair.
"""
function phycosphere_hein(config)
    produce_or_load(
        datadir("HeinMod"), config;
        prefix = "phycosphere",
        suffix = "jld2",
        tag = false,
        loadfile = false,
        force = true,
    ) do config
        @unpack R, Cₛ = jldopen(datadir("HeinMod", "RC.jld2"), "r")
        phycosphere_hein(config, R, Cₛ)
    end
end

function phycosphere_hein(config, R, Cₛ)
    @unpack C₀, T, Dc, U, Π = config
    S = map(p -> HeinModRadius(p..., C₀, T, Dc, U, Π), Iterators.product(R, Cₛ))
    @strdict S
end

## parameters
T = [50, 100, 200]u"ms"
Dc = [100, 500, 1000]u"μm^2/s"
U = [25, 50, 100]u"μm/s"
Π = [1, 6]u"1"
C₀ = [1]u"nM"
allparams = @strdict T Dc U Π C₀
dicts = dict_list(allparams)

## run
map(phycosphere_hein, dicts)
