using Test, LinearAlgebra, Random

include(joinpath(@__DIR__, "..", "src", "PaperAllToAll.jl"))
using .PaperAllToAll

module ScalingLegacy
include(joinpath(@__DIR__, "..", "LRSD Scaling Analysis"))
end

module BellLegacy
include(joinpath(@__DIR__, "..", "Stabilizer Nullity (Magic) Computation"))
end

@testset "All-to-all dense reference" begin
    cliffords = two_qubit_cliffords()
    @test length(cliffords) == 11520
    @test all(norm(U' * U - I) < 1e-10 for U in cliffords[1:100:end])

    for pt in (0.0, 0.3)
        trace = purification_trajectory(3; p_m=0.0, p_t=pt, tf=6,
                                         rng=MersenneTwister(9), cliffords)
        @test all(isapprox(s, log(2); atol=1e-10) for s in trace)
    end

    measured = ComplexF64[0, 0, 0, 1]
    event = paper_step!(MersenneTwister(1), measured, 2;
                        p_cz=0.0, p_t=0.0, p_m=1.0, basis=:Z)
    bit = 1 << (event.site - 1)
    p0 = sum(abs2(measured[i + 1]) for i in 0:3 if (i & bit) == 0)
    x_expect = sum(2real(conj(measured[i + 1]) * measured[(i | bit) + 1])
                   for i in 0:3 if (i & bit) == 0)
    @test event.outcome == 1
    @test isapprox(p0, 0.5; atol=1e-12)
    @test isapprox(x_expect, 1.0; atol=1e-12)

    fast_purification = purification_ensemble(3; p_m=1.0, tf=18,
                                               trajectories=8, seed=1, cliffords)
    @test fast_purification.mean[end] < 1e-10

    ensemble = purification_ensemble(3; p_m=0.4, tf=6, trajectories=4,
                                      seed=12, cliffords)
    repeated = purification_ensemble(3; p_m=0.4, tf=6, trajectories=4,
                                      seed=12, cliffords)
    @test ensemble.samples == repeated.samples
    @test size(ensemble.samples) == (7, 4)
    @test all(0 <= s <= log(2) + 1e-10 for s in ensemble.mean)

    t = collect(0:30)
    y = exp.(-t ./ 7)
    fit = fit_late_tail(t, y, zeros(length(t)), 4)
    @test fit !== nothing
    @test isapprox(fit.tau, 7.0; atol=1e-10)
    @test isapprox(fit.r2, 1.0; atol=1e-10)
end

@testset "Original LRSD and Bell kernels" begin
    for px in (0.0, 1.0)
        Random.seed!(123)
        result = ScalingLegacy.compute_logno(2, 0.8, px, 0.2, 1.0,
                                              0.5, 1e-15, 4, false, "")
        @test length(result.logical_terms) == 2
        @test all(1 <= k <= 256 for k in result.logical_terms)
    end

    Random.seed!(23)
    @test BellLegacy.runy(4, 0.0, 0.5, 1.0, 2, "") == [0.0, 0.0]
    Random.seed!(23)
    magic = BellLegacy.runy(4, 1.0, 0.5, 1.0, 2, "")
    @test length(magic) == 2
    @test all(0 <= m <= 4 for m in magic)
end
