include("apply_noise.jl")
include("apply_dyson.jl")

    """
        (p, t, pmax, tmax) = spatialsearch(psi0, Adjacency; kwargs...)

    Evaluate a (noisy) spatial search run on the graph described by the adjacency
    matrix Adjacency, starting from an initial state `psi0`. Return a `NamedTuple`
    consisting of
    - `p`: the probability of finding the walker in the target node as a
    f unction of `t`
    - `t`: the time instants at which probArray is evaluated
    - `pmax`: the maximum probability
    - `tmax`: the optimal time at which `pmax` is reached for the first time

    This function implements the spatial search as devised by Childs & Goldstone,
    Phys. Rev. A 70, 22314 (2004), with the addition of RTN noise as described in
    the manuscript

    # Arguments
    ## Required
    - `psi0`: a complex vector representing the initial state
    - `Adjacency`: a function that returns a sparse adjacency matrix for the desired
        graph topology

    ## Keyword arguments
    - `maxdt::Real = 0.05`: maximum allowed dt for the output
    - `gamma::Real = 1.`: Coupling γ in the Hamiltonian
    - `posW::Integer = 1`: Position of the target node
    - `time::Real = 1.`: Final time
    - `mu::Real = 1.`: Switching rate μ of the noise
    - `noiseStrength::Real = 1.`: Coupling ν to the noise (noise stregth)
    - `noiseRealizations::Integer = 0`: Number of Montecarlo noise realizations
    - `dysonOrder::Integer = 4`: Order of the Dyson series expansion

    # Examples

    In this example we run the noiseless spatial search on a complete graph with
    `N = 4` and we check that `p_max = 1` and `t_opt = 1` (up to error due to
    `dt`)

    ```jldoctest
    julia> (probArray, t, p_max, t_opt) = spatialsearch(superposition_state(4),
    complete_graph_Ad; time=5, gamma = 1. / 4, maxdt = .001)

    julia> isapprox(p_max, 1, atol=1e-6)
    true

    julia> isapprox(t_opt, pi, atol=1e-3)
    true
    ```
    """
    function spatialsearch(psi0, Adjacency;
                        maxdt::Real=0.05, # maximum allowed time for dt
                        gamma::Real=1.,
                        posW::Integer=1,
                        time::Real=1.,
                        mu::Real=1.,
                        coupling::Real=1.,
                        noiseStrength::Real=0.,
                        noiseRealizations::Integer=1,
                        dysonOrder::Integer=4)

        # maxdt = maximum value of dt of the evolution.
        # posW = position of the state we want to find in the basis of the nodes
        # of the graph (position in the state array).

        N = length(psi0)   # N (dimension of the lattice)

        # Required number of output timesteps
        output_timesteps::Int = Int(cld(time, maxdt) + 1)

        # We evaluate the number of internal steps required so that the noise is
        # sampled correctly. kout is the number of output timesteps t
        (Mtot, kout) = internal_M(time, mu, output_timesteps)

        # every internal step has a length = dt
        dtev = internal_dt(time, mu, output_timesteps)

        Ad = Adjacency(N)  # We define the whole adjacency matrix

        # We obtain the Hamiltonian (Laplacian matrix) from the adjacency matrix
        H = discLapl(Ad)

        # The number of links is the number of nonzero elements in the upper
        # triangular adjacency matrix.
        link_number = Int(length(Ad.nzval) / 2)
        tvec = range(0., stop=time, length=Mtot)

        # Contains the output
        probArray = zeros(Mtot)
        #probArray = @distributed (+) for n = 1 : noiseRealizations
        for n = 1 : noiseRealizations
            # array containing the probability of getting the required state |posW>
            # if we measure the evoluted state, at each internal timestep t
            tmp = zeros(Float64, Mtot)
            tmp[1] = abs2(psi0[posW]) / noiseRealizations

            psi = copy(psi0)  #psi0 is the initial state (array of complex floats)
            #kq = similar(psi)
            # This is the noise array generated for each indipedent link.
            # It's a Mtot x NumberOfIndipendentLinks matrix.
            temparray = - gamma * (coupling .+
                        noiseStrength .* generateRTN(tvec, mu, link_number))

            for t = 2 : Mtot    # evolution over time
                apply_noise!(H, Ad, view(temparray, t-1, :))

                # we subtract the oracle Hamiltonian from the noisy Laplacian
                H[posW, posW] -= 1.

                #copyto!(kq, psi)

                # evolution (dysonOrder is the maximum order of dt)
                # for k = 1 : dysonOrder
                #     kq = - 1im * dtev/k * H * kq
                #     psi += kq
                # end
                apply_dyson!(H, psi, dtev, dysonOrder)

                psi ./= norm(psi)     # state psi at time=t

                # we calculate the probability
                tmp[t] = abs2(psi[posW]) ./ noiseRealizations

            end  # end of loop over time

            probArray += tmp # The parallel-for reduction sums this value to probArray
        end  # end of loop over noise realizations

        # we find the maximum and its position in the probability array
        (max, position) = findmax(probArray)

        return  (p=probArray, t=tvec, pmax=max, tmax=(position - 1) * dtev)
    end
