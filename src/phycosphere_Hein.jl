export HeinRadius

"""
    HeinRadius(R, L, C₀, T, U, Π; q=1)
Find the distance from the source at which the SNR drops below the threshold.

*Parameters*
- `R`: source radius
- `L`: source leakage rate
- `C₀`: background concentration
- `T`: sensory integration timescale
- `U`: convective speed of the sensor
- `Π`: chemotactic precision factor
- `q`: SNR threshold, defaults to 1
"""
function HeinRadius(R, L, C₀, T, U, Π; q=1)
    # equation to solve: SNR = q
    f(r) = snr(r,L,C₀,T,U,Π) - q
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