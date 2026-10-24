#!/usr/bin/env julia

using CSV, DataFrames, Random
include(joinpath(@__DIR__, "cli.jl"))
using .ExampleCLI
include(joinpath(@__DIR__, "..", "Stabilizer Nullity (Magic) Computation"))

opts = options(ARGS, Dict(
    "L" => "4,6", "pm" => "0.5", "eta" => "1", "beta" => "1",
    "trajectories" => "4", "seed" => "2026",
    "out" => joinpath(@__DIR__, "..", "outputs", "fig6_bell_magic.csv")))

sizes, pms = int_values(opts["L"]), float_values(opts["pm"])
eta, betas = parse(Float64, opts["eta"]), float_values(opts["beta"])
trajectories, seed = parse(Int, opts["trajectories"]), parse(Int, opts["seed"])
all(L -> L >= 2, sizes) || throw(ArgumentError("all sizes must be >= 2"))
all(pm -> 0 <= pm <= 1, pms) || throw(ArgumentError("pm must lie in [0,1]"))
trajectories > 0 || throw(ArgumentError("trajectories must be positive"))

rows = NamedTuple[]
for (li, L) in enumerate(sizes), (pi, pm) in enumerate(pms), (bi, beta) in enumerate(betas)
    pt = eta / L^beta
    0 <= pt <= 1 || throw(ArgumentError("eta/L^beta must lie in [0,1]"))
    point_seed = seed + 100000 * li + 10000 * pi + 1000 * bi
    Random.seed!(point_seed)
    magic = runy(L, eta / L, pm, beta, trajectories, "")
    for trial in 1:trajectories
        push!(rows, (L=L, pm=pm, eta=eta, beta=beta, pt=pt, seed=point_seed,
                     trajectory=trial, bell_samples=4L, M_estimate=magic[trial]))
    end
end

out = abspath(opts["out"])
mkpath(dirname(out))
CSV.write(out, DataFrame(rows))
println("wrote $(length(rows)) trajectory rows to $out")
