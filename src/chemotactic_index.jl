export IC

"""
    Kn(λ,R) = 4λ/3R
Knudsen number for a searcher with persistence length `λ`
and a target with radius `R`.
"""
@inline Kn(λ,R) = 4λ/3R

"""
    IC(λ,R,S)
Evaluate the black-hole limit of the chemotactic index for a searcher with
persistence length `λ`, a target of radius `R` and a phycosphere of radius `S`.
"""
@inline IC(λ,R,S) = Kn(λ,R)/Kn(λ,S) * (1+Kn(λ,R))/(1+Kn(λ,S))