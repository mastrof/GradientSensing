export nevents, eventrate, waitingtime

"""
    nevents(c, T, Dc)
Evaluate average number of counting events occuring during a time interval `T`
when the local concentration is `c` and the diffusivity of the compound is `Dc`.
Uses the global defaults for the sensor radius `a=0.5μm`. 
"""
@inline nevents(c, T, Dc) = upreferred(4π*Dc*a*c*T*Unitful.Na)

"""
    eventrate(c)
Evaluate average rate of counting events when the local concentration is `c`
and the diffusivity of molecules is `Dc`.
Uses the global defaults for the sensor radius `a=0.5μm`. 
"""
@inline eventrate(c, Dc) = upreferred(4π*Dc*a*c*Unitful.Na)

"""
    waitingtime(c)
Evaluate average waiting time between events when the local concentration is `c`
and the molecular diffusivity is `Dc`.
Uses the global defaults for the sensor radius `a=0.5μm`. 
"""
@inline waitingtime(c, Dc) = 1/eventrate(c, Dc)
