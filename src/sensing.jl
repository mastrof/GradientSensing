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
@inline function σMW(r, L, C₀, T)
    Cr = C(r, L, C₀)
    sqrt(3*Cr / (π*a*Dc*T^3*Unitful.Na)) |> u"μM/s"
end

"""
    signal(r, L, U)
True signal from a temporal measurement of a spatial gradient
(equivalent to a convective derivative).
"""
@inline signal(r, L, U) = U*∇C(r,L) |> u"μM/s"

"""
    noise(r, L, C₀, T, Π)
Noise on gradient sensing obtained multiplying the Mora-Wingreen noise `σMW`
by the chemotactic precision factor `Π`.
"""
@inline noise(r, L, C₀, T, Π) = Π * σMW(r, L, C₀, T) |> u"μM/s"

"""
    snr(r, L, C₀, T, U, Π)
Signal to noise ratio.

*Parameters*
- `r`: distance from the center of the source
- `L`: leakage rate of the source
- `C₀`: background concentration (at infinite distance from the source)
- `T`: sensory integration timescale
- `U`: convective speed of the sensor
- `Π`: chemotactic precision factor
"""
@inline snr(r, L, C₀, T, U, Π) = signal(r,L,U)/noise(r,L,C₀,T,Π)