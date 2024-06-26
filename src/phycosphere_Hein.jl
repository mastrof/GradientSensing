export HeinRadius, HeinModRadius, HeinFDRadius, HeinRFRadius

"""
    HeinRadius(R, L, C₀, T, Dc, U, Π; q=1)
Find the distance from the source at which the SNR drops below the threshold.

*Parameters*
- `R`: source radius
- `L`: source leakage rate
- `C₀`: background concentration
- `T`: sensory integration timescale
- `Dc`: diffusivity of compound
- `U`: convective speed of the sensor
- `Π`: chemotactic precision factor
- `q`: SNR threshold, defaults to 1
"""
function HeinRadius(
    R::Quantity{<:Real,𝐋},
    L::Quantity{<:Real,𝐍 / 𝐓},
    C₀::Quantity{<:Real,𝐍 / 𝐋^3},
    T::Quantity{<:Real,𝐓},
    Dc::Quantity{<:Real,𝐋^2 / 𝐓},
    U::Quantity{<:Real,𝐋 / 𝐓},
    Π;
    q=1
)
    # equation to solve: SNR = q
    f(r) = snr(r, L, C₀, T, Dc, U, Π) - q
    # try to solve SNR=q
    # if a solution cannot be found, there is no phycosphere
    # so h is set to h=R
    h = try
        find_zero(f, R)
    catch e
        R
    end
    # if the estimated phycosphere is smaller than R, set it to R
    h > R ? h : R
end

"""
    HeinRadius(R, Cₛ, C₀, T, Dc, U, Π; q=1)
Find the distance from the source at which the SNR drops below the threshold.

*Parameters*
- `R`: source radius
- `Cₛ`: concentration at the surface of the source
- `C₀`: background concentration
- `T`: sensory integration timescale
- `Dc`: diffusivity of compound
- `U`: convective speed of the sensor
- `Π`: chemotactic precision factor
- `q`: SNR threshold, defaults to 1
"""
function HeinRadius(
    R::Quantity{<:Real,𝐋},
    Cₛ::Quantity{<:Real,𝐍 / 𝐋^3},
    C₀::Quantity{<:Real,𝐍 / 𝐋^3},
    T::Quantity{<:Real,𝐓},
    Dc::Quantity{<:Real,𝐋^2 / 𝐓},
    U::Quantity{<:Real,𝐋 / 𝐓},
    Π;
    q=1
)
    # equation to solve: SNR = q
    f(r) = snr(r, R, Cₛ, C₀, T, Dc, U, Π) - q
    # try to solve SNR=q
    # if a solution cannot be found, there is no phycosphere
    # so h is set to h=R
    h = try
        find_zero(f, R)
    catch e
        R
    end
    # if the estimated phycosphere is smaller than R, set it to R
    h > R ? h : R
end

"""
    HeinModRadius(R, Cₛ, C₀, T, Dc, U, Π; q=1)
Find the distance from the source at which the modulated SNR drops below the threshold.
"""
function HeinModRadius(R, Cₛ::Quantity{<:Real,𝐍 / 𝐋^3}, C₀, T, Dc, U, Π; q=1)
    λ = U * T
    f(r) = snr(r, R, Cₛ, C₀, T, Dc, U, Π) * modulator((r / λ)^1.5) - q
    #f(r) = snr(r, R, Cₛ, C₀, T, Dc, U, Π) * fp(r / λ) - q
    h = try
        find_zero(f, R)
    catch e
        R
    end
    h > R ? h : R
end

@inline modulator(x) = 1 - exp(-x)
@inline fp(x) = min((x / 2)^1, 1.0)
@inline rf(x) = sqrt(1 + 3/20 * x^2)

"""
    HeinFDRadius(R, Cₛ, C₀, T, Dc, U, Π; q=1)
Find the distance from the source at which the SNR, with the signal estimated
via finite differences, drops below the threshold.
"""
function HeinFDRadius(R, Cₛ::Quantity{<:Real,𝐍 / 𝐋^3}, C₀, T, Dc, U, Π; q=1)
    λ = U * T
    f(r) = snrFD(r, R, Cₛ, C₀, T, Dc, U, Π) - q
    #f(r) = snr(r, R, Cₛ, C₀, T, Dc, U, Π) * (1 / (1 + λ / r)) - q
    h = try
        find_zero(f, 2R, Order16(); xatol=1e-16u"μm")
    catch e
        #println(e)
        R
    end
    h > R ? h : R
end

"""
    HeinRFRadius(R, Cₛ, C₀, T, Dc, U, Π; q=1)
Find the distance from the source at which the SNR, with the noise estimated
assuming a 2nd order expansion of the concentration field.
"""
function HeinRFRadius(R, Cₛ::Quantity{<:Real,𝐍/𝐋^3}, C₀, T, Dc, U, Π; q=1)
    λ = U * T
    f(r) = snr(r, R, Cₛ, C₀, T, Dc, U, Π) / rf(λ/r) - q
    h = try
        find_zero(f, R, Order0(); atol=1e-3u"μm")
    catch e
        R
    end
    h > R ? h : R
end
