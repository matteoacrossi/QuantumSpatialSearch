Profile.clear()

include("../Grover/main.jl")
include("../Grover/evolutionGrover.jl")
include("../Grover/RTN.jl")
include("../Grover/graphs.jl")

N = 10
maxdt = 0.01
psi0 = 1/sqrt(N) * ones(ComplexF64, N)

@profile spatialsearch(psi0, star_graph_Ad;
                    maxdt=maxdt, # maximum allowed time for dt
                    gamma=1./N,
                    posW=1,
                    time=2*sqrt(N),
                    mu=1.,
                    coupling=2.,
                    noiseStrength=2.,
                    noiseRealizations=100,
                    dysonOrder=4);
