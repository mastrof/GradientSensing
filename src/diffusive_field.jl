export C, ∇C

"""
    C(r, L, C₀)
Concentration at distance `r` from the center of a spherical diffusive source
with leakage rate `L` and background concentration `C₀`:
C(r) = C₀ + L/(4π*Dc*r).

Returns the result in μM.
"""
@inline C(r, L::Quantity{<:Real,𝐍/𝐓}, C₀) = C₀ + L/(4π*Dc*r)  |> u"μM"
"""
    C(r, R, Cₛ, C₀)
Concentration at distance `r` from the center of a spherical diffusive source
of radius `R` with surface concentration `Cₛ` and background concentration `C₀`:
C(r) = C₀ + Cₛ*R/r.

Returns the result in μM.
"""
@inline C(r, R, Cₛ::Quantity{<:Real,𝐍/𝐋^3}, C₀) = C₀ + Cₛ*R/r |> u"μM"
"""
    ∇C(r, L)
Absolute value of the concentration gradient at distance `r` from the center
of a spherical diffusive source with leakage rate `L`:
∇C(r) = L/(4π*Dc*r^2).

Returns the result in "μM/μm".
"""
@inline ∇C(r, L::Quantity{<:Real,𝐍/𝐓}) = L/(4π*Dc*r^2) |> u"μM/μm"
"""
    ∇C(r, R, Cₛ)
Absolute value of the concentration gradient at distance `r` from the center
of a spherical diffusive source of radius `R` with surface concentration `Cₛ`:
∇C(r) = Cₛ*R/r^2.

Returns the result in "μM/μm".
"""
@inline ∇C(r, R, Cₛ::Quantity{<:Real,𝐍/𝐋^3}) = Cₛ*R/r^2 |> u"μM/μm"