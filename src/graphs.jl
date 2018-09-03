using SparseArrays

"""
    complete_graph_Ad(N)

Returns a sparse adjacency matrix for the complete graph with `N` nodes.

```jldoctest
julia> Matrix(complete_graph_Ad(N))
3×3 Array{Float64,2}:
 0.0  1.0  1.0
 1.0  0.0  1.0
 1.0  1.0  0.0
```
"""
function complete_graph_Ad(N)
    return dropzeros(1. .- sparse(1.0I, N, N))
end;

"""
    star_graph_Ad(N)

Returns a sparse adjacency matrix for the complete graph with `N` nodes.

```jldoctest
julia> Matrix(star_graph_Ad(3))
3×3 Array{Float64,2}:
 0.0  1.0  1.0
 1.0  0.0  0.0
 1.0  0.0  0.0
```
"""
function star_graph_Ad(N)
    V = Array{Int64}(undef, 2N-2)
    U = Array{Int64}(undef, 2N-2)
    L = Array{Float64}(undef, 2N-2)

    for i=1:N-1
        V[i]=1
        U[i]=i+1
        L[i]=1.
    end

    for i=1:N-1
        V[i+N-1]=i+1
        U[i+N-1]=1
        L[i+N-1]=1.
    end

    return sparse(V,U,L)
end

"""
    discLapl(adj::SparseMatrixCSC)

Return the Laplacian matrix corresponding to the adjacency matrix `adj`.

```jldoctest
julia> Matrix(discLapl(star_graph_Ad(3)))
3×3 Array{Float64,2}:
  2.0  -1.0  -1.0
 -1.0   1.0   0.0
 -1.0   0.0   1.0
```
"""
function discLapl(adj::SparseMatrixCSC)

    temp = spzeros(maximum(adj.rowval), maximum(adj.rowval))

    for i=1:maximum(adj.rowval)
        temp[i,i]=adj.colptr[i+1]-adj.colptr[i]
        for j=adj.colptr[i]:(adj.colptr[i+1]-1)
            temp[adj.rowval[j],i]=-adj[adj.rowval[j],i]
        end
    end

    return temp
end
