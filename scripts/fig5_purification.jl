#!/usr/bin/env julia

using CSV, DataFrames
include(joinpath(@__DIR__, "cli.jl"))
using .ExampleCLI
include(joinpath(@__DIR__, "..", "src", "PaperAllToAll.jl"))
using .PaperAllToAll

opts = options(ARGS, Dict(
    "L" => "4", "pm" => "0.2,0.4", "eta" => "0", "beta" => "2",
    "trajectories" => "16", "seed" => "2026", "tf" => "auto",
    "out" => joinpath(@__DIR__, "..", "outputs", "fig5_purification.csv")))

sizes, pms = int_values(opts["L"]), float_values(opts["pm"])
eta, beta = parse(Float64, opts["eta"]), parse(Float64, opts["beta"])
trajectories, seed = parse(Int, opts["trajectories"]), parse(Int, opts["seed"])
all(L -> 2 <= L <= 12, sizes) || throw(ArgumentError("dense reference needs 2 <= L <= 12"))
all(pm -> 0 <= pm <= 1, pms) || throw(ArgumentError("pm must lie in [0,1]"))
trajectories > 0 || throw(ArgumentError("trajectories must be positive"))

cliffords = two_qubit_cliffords()
rows = NamedTuple[]
for (li, L) in enumerate(sizes), (pi, pm) in enumerate(pms)
    pt = eta / L^beta
    0 <= pt <= 1 || throw(ArgumentError("eta/L^beta must lie in [0,1]"))
    tf = opts["tf"] == "auto" ? 2L^2 : parse(Int, opts["tf"])
    point_seed = seed + 10000 * li + 1000 * pi
    result = purification_ensemble(L; p_m=pm, p_t=pt, tf, trajectories,
                                    seed=point_seed, cliffords)
    fit = fit_late_tail(result.time, result.mean, result.sem, L)
    for t in eachindex(result.time)
        push!(rows, (L=L, pm=pm, eta=eta, beta=beta, pt=pt,
                     trajectories=trajectories, seed=point_seed, t=result.time[t],
                     mean_SQ=result.mean[t], sem_SQ=result.sem[t],
                     tau=fit === nothing ? missing : fit.tau,
                     fit_start=fit === nothing ? missing : fit.t_start,
                     fit_end=fit === nothing ? missing : fit.t_end,
                     fit_R2=fit === nothing ? missing : fit.r2))
    end
end

out = abspath(opts["out"])
mkpath(dirname(out))
CSV.write(out, DataFrame(rows))
println("wrote $(length(rows)) time-series rows to $out")
