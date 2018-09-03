module QuantumSpatialSearch

    using SparseArrays
    using LinearAlgebra

    export spatialsearch

    const jumpProb = 0.05;   #jump probability in a timestep

    include("graphs.jl")
    include("spatialsearch.jl")
    include("utilities.jl")
    include("RTN.jl")

end # module
