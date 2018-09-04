
"""
    internal_M(time, mu, timesteps)

Return `(numberOfSteps, kout)` where `numberOfSteps` is the required number of
internal simulation steps so that the probability of having a noise jump in a
timestep is `jumpProb`. `kout` tells how many internal time steps are there
between timesteps with output.

`timesteps` = number of output steps.
"""
function internal_M(time, mu, timesteps)
    temp=convert(Int64,cld(time,jumpProb/mu )+1) # minimum number of internal steps

    if time/(temp-1)>0.1
        temp=ceil(Int,time/0.1+1)
    end
    if temp<=timesteps                            # if temp<timesteps then I simply use the number of output steps
        kout=1
        return timesteps, kout
    else

    # otherwise we choose the smallest multiple of (timesteps-1) bigger than the required number of internal steps -1 :
    # numberOfSteps=kout*(timesteps-1)+1

        kout=cld((temp-1),(timesteps-1))
        # in finding kout we have to use (temp-1) and (timesteps-1) because the time step=1 is fixed for both
        # internal and output steps.

        return ((timesteps-1)*kout+1),kout

    end
end


internal_dt(time,mu,timesteps)=time/(internal_M(time,mu,timesteps)[1]-1) # unit of time given the number of internal steps


rescaledMu(time,mu,timesteps)=dt(time,mu,timesteps)*mu


"""
    psi = superposition_state(N::Integer)

Return the superposition state for the `N`-dimensional Hilbert space. `psi` is
a `Array{Complex{Float64}}`.

```jldoctest
julia> superposition_state(4)
4-element Array{Complex{Float64},1}:
 0.5 + 0.0im
 0.5 + 0.0im
 0.5 + 0.0im
 0.5 + 0.0im
```
"""
function superposition_state(N::Integer)
    return ones(ComplexF64, N) ./ sqrt(N)
end
