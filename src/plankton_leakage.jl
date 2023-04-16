export leakage_rate

@inline amount_of_substance(R::Quantity) = bⱼ/nC * (R |> u"cm")^2.28
"""
    leakage_rate(R, PER)
Evaluate the leakage rate from a spherical phytoplankton cell of radius `R`
and percent extracellular release `PER`.
"""
@inline leakage_rate(R, PER) = μ*PER*amount_of_substance(R) |> upreferred
