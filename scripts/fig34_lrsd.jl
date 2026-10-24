#!/usr/bin/env julia

using CSV, DataFrames, Random
include(joinpath(@__DIR__, "cli.jl"))
using .ExampleCLI
include(joinpath(@__DIR__, "..", "LRSD Scaling Analysis"))

opts = options(ARGS, Dict(
    "L" => "4,6", "pm" => "0.6,0.8", "basis" => "X", "eta" => "1",
    "beta" => "1", "cutoff" => "1e-15", "trajectories" => "4",
    "seed" => "2026", "out" => joinpath(@__DIR__, "..", "outputs", "fig34_lrsd.csv")))

sizes = int_values(opts["L"])
pms = float_values(opts["pm"])
cutoffs = float_values(opts["cutoff"])
basis = uppercase(opts["basis"])
basis in ("X", "Z") || throw(ArgumentError("basis must be X or Z"))
etas, beta = float_values(opts["eta"]), parse(Float64, opts["beta"])
trajectories, seed = parse(Int, opts["trajectories"]), parse(Int, opts["seed"])
all(L -> L >= 2, sizes) || throw(ArgumentError("all sizes must be >= 2"))
all(pm -> 0 <= pm <= 1, pms) || throw(ArgumentError("pm must lie in [0,1]"))
all(c -> c >= 0, cutoffs) || throw(ArgumentError("cutoff must be nonnegative"))
trajectories > 0 || throw(ArgumentError("trajectories must be positive"))

rows = NamedTuple[]
for (li, L) in enumerate(sizes), (pi, pm) in enumerate(pms), (ei, eta) in enumerate(etas)
    pt = eta / L^beta
    0 <= pt <= 1 || throw(ArgumentError("eta/L^beta must lie in [0,1]"))
    paired_seed = seed + 100000 * li + 10000 * pi + 1000 * ei
    for cutoff in cutoffs
        Random.seed!(paired_seed)
        result = compute_logno(trajectories, pm, basis == "X" ? 1.0 : 0.0,
                               eta / L, beta, 0.5, cutoff, L, false, "")
        for trial in 1:trajectories
            k = Int(result.logical_terms[trial])
            push!(rows, (L=L, pm=pm, basis=basis, eta=eta, beta=beta, pt=pt,
                         cutoff=cutoff, seed=paired_seed, trajectory=trial,
                         logical_terms=k, representation_entries=(2L + 1)^2 + (L + 1) * k,
                         seconds=result.seconds[trial]))
        end
    end
end

out = abspath(opts["out"])
mkpath(dirname(out))
CSV.write(out, DataFrame(rows))
println("wrote $(length(rows)) trajectory rows to $out")
