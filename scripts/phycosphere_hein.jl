using DrWatson
@quickactivate :GradientSensing
using JLD2

function phycosphere_hein()
    f = jldopen(datadir("Hein", "paramspaceRL.jld2"))
    R, L = f["R"], f["L"]
    
    # set C₀, T, U, Π to some default values
    # can't assign units now or they are ignored in the output filename
    C₀ = 1#u"nM"
    T = 100#u"ms"
    U = 46.5#u"μm/s"
    Π = 6
    params = @strdict C₀ T U Π
    
    produce_or_load(
        datadir("Hein"), params;
        prefix="phycosphere", suffix="jld2",
        tag = false
    ) do params
        # reassign units to parameters
        C₀ = params["C₀"]u"nM"
        T = params["T"]u"ms"
        U = params["U"]u"μm/s"
        Π = params["Π"]u"1"
        # evaluate phycosphere radius for each (R,L) pair
        S = map(p -> HeinRadius(p..., C₀, T, U, Π), Iterators.product(R,L))
        @strdict S
    end
end

phycosphere_hein()