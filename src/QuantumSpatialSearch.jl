module QuantumSpatialSearch

    using SparseArrays
    using LinearAlgebra
    using Distributed
    
    export spatialsearch, generateRTN

    const jumpProb = 0.05;   #jump probability in a timestep

    include("graphs.jl")
    include("spatialsearch.jl")
    include("utilities.jl")
    include("RTN.jl")

end # module
