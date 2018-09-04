using Distributed
using Test

@testset "Noiseless" begin
# In this test we check that the noiseless case gets to probability 1 in the
# expected time (up to tolerance depending on maxdt)
    N = 20
    maxdt = 0.01
    psi0 = 1. / sqrt(N) * ones(ComplexF64, N)
    @testset "Complete graph" begin
        res = spatialsearch(psi0, complete_graph_Ad;
                    maxdt=maxdt, # maximum allowed time for dt
                    gamma=1. / N,
                    posW=1,
                    time=2*sqrt(N),
                    mu=1.,
                    coupling=1.,
                    noiseStrength=0.,
                    noiseRealizations=1,
                    dysonOrder=4);

        @test res[1] >= zero(res[1]) # Positive probability
        @test res[3] ≈ 1 atol = 1e-5 # Maximum is 1
        @test res[4] ≈ pi*sqrt(N)/2  atol = maxdt # For the correct time

    end

    @testset "Star graph" begin
        @testset "Central target" begin
            res = spatialsearch(psi0, star_graph_Ad;
                        maxdt=maxdt, # maximum allowed time for dt
                        gamma=1. / N,
                        posW=1,
                        time=2*sqrt(N),
                        mu=1.,
                        coupling=1.,
                        noiseStrength=0.,
                        noiseRealizations=1,
                        dysonOrder=4);

            @test res[1] >= zero(res[1]) # Positive probability
            @test res[3] ≈ 1 atol = 1e-5
            @test res[4] ≈ pi*sqrt(N)/2 atol = maxdt
        end

        @testset "External target" begin
            res = spatialsearch(psi0, star_graph_Ad;
                        maxdt=maxdt, # maximum allowed time for dt
                        gamma=1.,
                        posW=2,
                        time=2*sqrt(N),
                        mu=1.,
                        coupling=1.,
                        noiseStrength=0.,
                        noiseRealizations=1,
                        dysonOrder=4);

            @test res[1] >= zero(res[1]) # Positive probability
            @test res[3] ≈ 1  atol = 1 / N^2
            @test res[4] ≈ pi*sqrt(N)/2  atol = maxdt
        end
    end
end

@testset "Noisy" begin
    N = 10
    maxdt = 0.01
    psi0 = 1/sqrt(N) * ones(ComplexF64, N)
    @testset "Star graph" begin
        @testset "Central target" begin
            res = spatialsearch(psi0, star_graph_Ad;
                        maxdt=maxdt, # maximum allowed time for dt
                        gamma=1. / N,
                        posW=1,
                        time=2*sqrt(N),
                        mu=0.01,
                        coupling=1.,
                        noiseStrength=1.,
                        noiseRealizations=1000,
                        dysonOrder=4);

            @test res[1] >= zero(res[1]) # Positive probability
            @test res[3] ≈ 0.5 atol = 5e-2
        end
    end
end


@testset "Distributed" begin
    addprocs(2)
    @everywhere using QuantumSpatialSearch
    N = 10
    maxdt = 0.01
    psi0 = 1/sqrt(N) * ones(ComplexF64, N)
    @testset "Star graph" begin
        @testset "Central target" begin
            res = spatialsearch(psi0, star_graph_Ad;
                        maxdt=maxdt, # maximum allowed time for dt
                        gamma=1. / N,
                        posW=1,
                        time=2*sqrt(N),
                        mu=0.01,
                        coupling=1.,
                        noiseStrength=1.,
                        noiseRealizations=1000,
                        dysonOrder=4);

            @test res[1] >= zero(res[1]) # Positive probability
            @test res[3] ≈ 0.5 atol = 5e-2
        end
    end
end
