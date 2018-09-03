#===========================================
Tests on the RTN generation routine function
===========================================#

#=
Here we calculate the autocorrelation functions of the two methods of generating RTN histories. The functions need a huge time in input (40000. is ok) in order to have ergodicity, and they return the autocorrelation array [<RTN(0)RTN(0)>-<RTN(0)><RTN(0)>,<RTN(1)RTN(0)>-<RTN(1)><RTN(0)>,...] up to the lag corresponding to an evolution time of tc=4./mu.
=#
using Test
using Statistics
using StatsBase

function autocorrelation(noisefunction, time, mu, dt)
    t = range(0., stop=time, step=dt)

    signal = noisefunction(t, mu)[:]
    tc = 4. / mu
    dt = Float64(t.step)

    lags = [x for x = 0:(Int(fld(tc, dt))-1)]
    autoC = autocov(signal, lags)

    tvec = range(0., stop=tc, length=Int(fld(tc,dt)))
    return (tvec, autoC)
end

function switchesAutoC(time, mu, timesteps)   # autocorrelation function of RTNgenSwitches

    time = range(0., stop=time, length=timesteps)
    temp = RTNgenSwitches(time, mu, 1)[:, 1]

    tc = 4. / mu
    dtt = Float64(time.step)
    lags = [x for x = 0:(Int(fld(tc, dtt))-1)]
    autoC = autocov(temp, lags)

    tvec = range(0., stop=tc, length=Int(fld(tc,dtt)))   # output "true" time from 0. to tc, split in the number of lags

    return (tvec,autoC)
end

function timeAutoC(time,mu,timesteps)       # autocorrelation function of RTNgenTime

    Mtot = M(time,mu,timesteps)[1]

    temp = RTNgenTime(time,mu,timesteps,1)[:,1]

    tc = 4. / mu
    dtt = time/Mtot   # time unit
    lags = [x for x=0:(Int(fld(tc,dtt))-1)]

    autoC = autocov(temp,lags)

    tvec = linspace(0., tc, Int(fld(tc,dtt)))   # output "true" time from 0. to tc, split in the number of lags

    return (tvec,autoC)
end

function spectrum(t, signal)
    n = length(signal)
    p = fft(signal)
    sampFreq = 1/(t[2] - t[1])
    nUniquePts = ceil(Int, (n+1)/2)
    p = p[1:nUniquePts]
    p = abs2.(p) / n^2
    # the fourier transform of the tone returned by the fft function contains both magnitude and phase information and is given in a complex representation (i.e. returns complex numbers). By taking the absolute value of the fourier transform we get the information about the magnitude of the frequency components.

    # odd nfft excludes Nyquist point
    if n % 2 > 0
        p[2 : length(p)] = p[2:length(p)]*2 # we've got odd number of   points fft
    else
        p[2: (length(p) -1)] = p[2: (length(p) -1)]*2 # we've got even number of points fft
    end

    freqArray = (0:(nUniquePts-1)) * (sampFreq / n)
    return (freqArray, p)
end


function test_autocorrelation(noisefunction, mu; atol=1e-1)
    (t, ac) = autocorrelation(noisefunction, 20000. /mu, mu, .05 / mu)
    @test ac ≈ exp.(-2 * mu * t) atol = atol
end

# Test that the mean is zero
@testset "RTN mean" begin
    @testset "γ = $γ" for γ in [0.01, 0.1 , 1., 10.]
        @test mean(generateRTN(0. : .05 / γ : 1000 / γ, γ, 100)) ≈ 0 atol=5e-3
    end
end

# Test that the autocorrelation function is the decaying exponential
@testset "RTN autocorrelation" begin
    @testset "γ = $γ" for γ in [0.01, 0.1 , 1., 10.]
        @testset "RTNgenTime" for i in 1:100
            test_autocorrelation(generateRTN, γ, atol=1e-1)
        end
    end
end
