    """
        (probArray, t, p_max, t_opt) = spatialsearch(psi0, Adjacency; kwargs...)

    Evaluate a (noisy) spatial search run on the graph described by the adjacency
    matrix Adjacency, starting from an initial state `psi0`. Return a tuple
    consisting of
    - `probArray`: the probability of finding the walker in the target node as a
    f unction of `t`
    - `t`: the time instants at which probArray is evaluated
    - `p_max`: the maximum probability
    - `t_opt`: the optimal time at which `p_max` is reached for the first time

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
    `N = 4` and we check that `p_max = 1` and `t_opt = 1` (up to error due to `dt`)

    ```jldoctest
    julia> (probArray, t, p_max, t_opt) = spatialsearch(superposition_state(4), complete_graph_Ad; time=5, gamma = 1. / 4, maxdt = .001)

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
        # posW = position of the state we want to find in the basis of the nodes of the graph (position in the state array). Usually we set         #        posW=1.

        N = length(psi0)   # N (dimension of the lattice)

        # Required number of output timesteps
        output_timesteps::Int = Int(cld(time, maxdt) + 1)

        # We evaluate the number of internal steps required so that the noise is
        # sampled correctly. kout is the number of output timesteps t
        (Mtot, kout) = M(time, mu, output_timesteps)

        # every internal step has a length = dt
        dtev = dt(time, mu, output_timesteps)

        Ad = Adjacency(N)  # We define the whole adjacency matrix
        H = discLapl(Ad)   # We obtain the Hamiltonian (discrete Laplacian) from the adjacency matrix

        # The number of links is the number of elements in the upper (lower) triangular adjacency matrix.
        link_number = Int(length(Ad.nzval) / 2)
        tvec = range(0., stop=time, length=Mtot)

        probArray = zeros(Mtot)

        for n = 1 : noiseRealizations

            # array containing the probability of getting the required state |posW>
            # if we measure the evoluted state, at each internal timestep t
            tmp = zeros(Float64, Mtot)
            tmp[1] = abs2(psi0[posW]) / noiseRealizations

            psi = copy(psi0)  #psi0 is the initial state (array of complex floats)

            # This is the noise array generated for each indipedent link.
            # It's a Mtot x NumberOfIndipendentLinks matrix.
            temparray = generateRTN(tvec, mu, link_number)

            for t = 2 : Mtot    # evolution over time

                # counter telling which of the indipendent links we are considering
                counter = 1

                for i = 1 : Int(sqrt(length(H)))

                    # this will be the value of the diagonal element [i,i],
                    # which depends on the noise on each
                    # link starting from the considered node i.
                    diagNoisyValue = 0.

                    # We perturb the laplacian matrix of the graph with the noise
                    # We do this at low level by acting on the sparse matrix
                    # representation

                    # elements of the row vector indicating which rows
                    # are occupied in column i
                    for j = Ad.colptr[i] : (Ad.colptr[i + 1] - 1)

                        # We insert "new" noise only in the lower triangular part
                        # of the Hamiltonian matrix
                        if i < Ad.rowval[j]

                            H[Ad.rowval[j], i] = - gamma * (coupling +
                                        noiseStrength .* temparray[t-1, counter])
                             # We have inserted noise in one of the indipedent links
                             counter+=1

                        elseif i > Ad.rowval[j]
                            H[Ad.rowval[j], i] = H[i, Ad.rowval[j]]
                        end

                        diagNoisyValue += -1. * H[Ad.rowval[j],i]
                    end

                    # the value of the diagonal element of a column (node i) must be = minus the sum of all
                    # the other values in the column, because of the conservation of the overall probability
                    # sum_j p(j|i)=1
                    H[i,i] = diagNoisyValue

                end

                H[posW, posW] -= 1.  # we subtract the oracle Hamiltonian from the noisy Laplacian

                kq = copy(psi)

                for k = 1 : dysonOrder      # evolution (dysonOrder is the maximum order of dt)
                    kq = - 1im * dtev/k * H * kq
                    psi += kq
                end

                psi ./= norm(psi)     # state psi at time=t

                # we calculate the probability
                tmp[t] = abs2(psi[posW]) ./ noiseRealizations

            end  # end of loop over time
            probArray += tmp # The parallel-for reduction sums this value to probArray
        end  # end of loop over noise realizations

        # we find the maximum and its position in the probability array
        (max, position) = findmax(probArray)

        # we return (in order) the probability array, the array of time, the maximum probability and the time at which we get the maximum
        return  (probArray, tvec, max, (position - 1) * dtev)
    end
