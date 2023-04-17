export nevents, eventrate, waitingtime

"""
    nevents(c, T)
Evaluate average number of counting events occuring during a time interval `T`
when the local concentration is `c`.
Uses the global defaults for the sensor radius `a=0.5μm` and
molecular diffusivity `Dc=608μm²/s`.
"""
@inline nevents(c, T) = upreferred(4π*Dc*a*c*T*Unitful.Na)

"""
    eventrate(c)
Evaluate average rate of counting events when the local concentration is `c`.
Uses the global defaults for the sensor radius `a=0.5μm` and
molecular diffusivity `Dc=608μm²/s`.
"""
@inline eventrate(c) = upreferred(4π*Dc*a*c*Unitful.Na)

"""
    waitingtime(c)
Evaluate average waiting time between events when the local concentration is `c`.
Uses the global defaults for the sensor radius `a=0.5μm` and
molecular diffusivity `Dc=608μm²/s`.
"""
@inline waitingtime(c) = 1/eventrate(c)