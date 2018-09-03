#Here we have two different methods of generating arrays of RTN.
using LinearAlgebra


function generateRTN(time::StepRangeLen, mu, noise_number::Int=1)
    # this one is the faster: we generate the time interval after which
    # we switch the RTN value
    timedim = time.len

    temp = Array{Float64}(undef, timedim, noise_number)
    rMu = mu * Float64(time.step) # Rescaled mu

    if rMu > jumpProb
        @warn "The given timestep is too large compared to the switching rate. Some jumps may be missed."
    end


    for n=1:noise_number    # we generate different values of RTN for each indipendent link of the graph

        i::Integer=1    # index of the noise array
        t::Integer=1    # starting time in units of timedim

        s = 1. * sign(1. - rand() * 2.)  # initial RTN value (either 1 or -1 with probability=0.5)
        # Notice that sign(0.) = 0 but here we can't get it since rand() is in [0, 1)

        # we push forward the time of a random value given by t=-floor(log(p0)/rescaledMu)
        # where p0 is a random number in [0,1], so that t is a random variable distributed
        # according to p(t)=Exp(-rescaledMu t) (no switches)

        t += Int(fld(-log(1. - rand()), rMu))

        while t < timedim

            view(temp, i:t, n) .= s      # we fill the noise array until t

            i=t+1
            s=-s                         # then we switch value
            t += 1 + Int(fld(-log(1. - rand()),rMu))   # again we push forward the time
        end

        view(temp, i:timedim, n) .= s        # last filling

    end

    return temp
end
