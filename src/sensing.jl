export snr

"""
    σMW(r, L, C₀, T)
Molecular counting noise on temporal gradient sensing for a perfect absorbing sphere
(Mora & Wingreen, Phys Rev Lett 2010).
Returns the result in "μM/s".

*Parameters*
- `r`: distance from the center of the source
- `L`: leakage rate of the source
- `C₀`: background concentration (at infinite distance from the source)
- `T`: sensory integration timescale
"""
@inline function σMW(r, L, C₀, T, Dc)
    Cr = C(r, L, C₀, Dc)
    sqrt(3 * Cr / (π * a * Dc * T^3 * Unitful.Na)) |> u"μM/s"
end
"""
    σMW(r, R, Cₛ, C₀, T)
Molecular counting noise on temporal gradient sensing for a perfect absorbing sphere
(Mora & Wingreen, Phys Rev Lett 2010).
Returns the result in "μM/s".

*Parameters*
- `r`: distance from the center of the source
- `R`: radius of the source
- `Cₛ`: concentration at surface of the source
- `C₀`: background concentration (at infinite distance from the source)
- `T`: sensory integration timescale
- `Dc`: diffusivity of compound
"""
@inline function σMW(r, R, Cₛ, C₀, T, Dc)
    Cr = C(r, R, Cₛ, C₀)
    sqrt(3 * Cr / (π * a * Dc * T^3 * Unitful.Na)) |> u"μM/s"
end

"""
    signal(r, L, U, Dc)
True signal from a temporal measurement of a spatial gradient
(equivalent to a convective derivative).
"""
@inline signal(r, L, U, Dc) = U * ∇C(r, L, Dc) |> u"μM/s"
"""
    signal(r, R, Cₛ, U)
True signal from a temporal measurement of a spatial gradient
(equivalent to a convective derivative).
"""
@inline signal(r, R, Cₛ, U) = U * ∇C(r, R, Cₛ) |> u"μM/s"

finitediff_b(f, x, δ) = (f(x) - f(x + δ)) / δ # backward
#finitediff_b(f, x, δ) = (+3f(x) - 4f(x + δ / 2) + f(x + δ)) / δ # backward 2nd order
#finitediff_b(f, x, δ) = (f(x - δ / 2) - f(x + δ / 2)) / δ # central
@inline function signalFD(r, R, Cₛ, C₀, U, T)
    f(r) = C(r, R, Cₛ, C₀)
    dC = finitediff_b(f, r, U * T)
    U * dC |> u"μM/s"
end

"""
    noise(r, L, C₀, T, Dc, Π)
Noise on gradient sensing obtained multiplying the Mora-Wingreen noise `σMW`
by the chemotactic precision factor `Π`.
"""
@inline noise(r, L, C₀, T, Dc, Π) = Π * σMW(r, L, C₀, T, Dc) |> u"μM/s"
"""
    noise(r, R, Cₛ, C₀, T, Dc, Π)
Noise on gradient sensing obtained multiplying the Mora-Wingreen noise `σMW`
by the chemotactic precision factor `Π`.
"""
@inline noise(r, R, Cₛ, C₀, T, Dc, Π) = Π * σMW(r, R, Cₛ, C₀, T, Dc) |> u"μM/s"

"""
    snr(r, L, C₀, T, Dc, U, Π)
Signal to noise ratio.

*Parameters*
- `r`: distance from the center of the source
- `L`: leakage rate of the source
- `C₀`: background concentration (at infinite distance from the source)
- `T`: sensory integration timescale
- `Dc`: diffusivity of compound
- `U`: convective speed of the sensor
- `Π`: chemotactic precision factor
"""
@inline snr(r, L, C₀, T, Dc, U, Π) = signal(r, L, U, Dc) / noise(r, L, C₀, T, Dc, Π)
"""
    snr(r, R, Cₛ, C₀, T, Dc, U, Π)
Signal to noise ratio.

*Parameters*
- `r`: distance from the center of the source
- `R`: radius of the source
- `Cₛ`: concentration at surface of the source
- `C₀`: background concentration (at infinite distance from the source)
- `T`: sensory integration timescale
- `Dc`: diffusivity of compound
- `U`: convective speed of the sensor
- `Π`: chemotactic precision factor
"""
@inline snr(r, R, Cₛ, C₀, T, Dc, U, Π) = signal(r, R, Cₛ, U) / noise(r, R, Cₛ, C₀, T, Dc, Π)

@inline snrFD(r, R, Cₛ, C₀, T, Dc, U, Π) = signalFD(r, R, Cₛ, C₀, U, T) / noise(r, R, Cₛ, C₀, T, Dc, Π)
