using QuantumSpatialSearch
using BenchmarkTools

N = 10
psi0 = 1/sqrt(N) * ones(ComplexF64, N)
maxdt = 0.01
println("Precompiling")
spatialsearch(psi0, star_graph_Ad;
                    maxdt=maxdt, # maximum allowed time for dt
                    gamma=1. / N,
                    posW=1,
                    time=2*sqrt(N),
                    mu=0.01,
                    coupling=1.,
                    noiseStrength=1.,
                    noiseRealizations=2,
                    dysonOrder=4);

println("Running benchmark")
res = @benchmark spatialsearch(psi0, star_graph_Ad;
                    maxdt=maxdt, # maximum allowed time for dt
                    gamma=1. / N,
                    posW=1,
                    time=2*sqrt(N),
                    mu=0.01,
                    coupling=1.,
                    noiseStrength=1.,
                    noiseRealizations=100,
                    dysonOrder=4)
