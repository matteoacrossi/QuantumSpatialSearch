function apply_noise!(H::SparseMatrixCSC, Ad::SparseMatrixCSC, noisearray)

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

            # We insert noise only in the lower triangular part
            # of the Hamiltonian matrix
            if i < Ad.rowval[j]

                H[Ad.rowval[j], i] = noisearray[counter]

                 counter+=1

            elseif i > Ad.rowval[j]
                H[Ad.rowval[j], i] = H[i, Ad.rowval[j]]
            end

            diagNoisyValue += -1. * H[Ad.rowval[j],i]
        end

        # the value of the diagonal element of a column (node i)
        # must be = minus the sum of all the other values in the
        # column, because of the conservation of the overall
        # probability sum_j p(j|i)=1
        H[i,i] = diagNoisyValue
    end
end
