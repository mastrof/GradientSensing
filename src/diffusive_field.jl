export C, ∇C

"""
    C(r, L, C₀)
Concentration at distance `r` from the center of a spherical diffusive source
with leakage rate `L` and background concentration `C₀`:
C(r) = C₀ + L/(4π*Dc*r).

Returns the result in μM.
"""
@inline C(r, L, C₀) = C₀ + L/(4π*Dc*r)  |> u"μM"
"""
    ∇C(r, L)
Absolute value of the concentration gradient at distance `r` from the center
of a spherical diffusive source with leakage rate `L`:
∇C(r) = L/(4π*Dc*r^2).

Returns the result in "μM/μm".
"""
@inline ∇C(r, L) = L/(4π*Dc*r^2) |> u"μM/μm"