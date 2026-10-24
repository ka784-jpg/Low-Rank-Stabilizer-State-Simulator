module PaperAllToAll

using LinearAlgebra
using Random
using Statistics

export two_qubit_cliffords, paper_step!, ancilla_entropy, purification_trajectory,
       purification_ensemble, fit_late_tail

const H = ComplexF64[1 1; 1 -1] / sqrt(2)
const X = ComplexF64[0 1; 1 0]
const S = ComplexF64[1 0; 0 im]
const CNOT = ComplexF64[1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
const SWAP = ComplexF64[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
const I2 = Matrix{ComplexF64}(I, 2, 2)

function canonical_clifford(U)
    first_nonzero = findfirst(z -> abs(z) > 1e-9, U)
    V = U / (U[first_nonzero] / abs(U[first_nonzero]))
    return V
end

function clifford_key(U)
    V = canonical_clifford(U)
    return Tuple(round(Int, component * 1e8) for z in V for component in (real(z), imag(z)))
end

function two_qubit_cliffords()
    generators = (kron(H, I2), kron(I2, H), kron(S, I2), kron(I2, S),
                  CNOT, SWAP * CNOT * SWAP)
    group = Matrix{ComplexF64}[Matrix{ComplexF64}(I, 4, 4)]
    seen = Set([clifford_key(group[1])])
    head = 1
    while head <= length(group)
        U = group[head]
        for G in generators
            V = canonical_clifford(G * U)
            key = clifford_key(V)
            if !(key in seen)
                push!(seen, key)
                push!(group, V)
            end
        end
        head += 1
    end
    length(group) == 11520 || error("two-qubit Clifford enumeration is incomplete")
    return group
end

function single_gate!(psi, U, q)
    bit = 1 << (q - 1)
    for i in 0:(length(psi) - 1)
        (i & bit) == 0 || continue
        a, b = psi[i + 1], psi[(i | bit) + 1]
        psi[i + 1] = U[1, 1] * a + U[1, 2] * b
        psi[(i | bit) + 1] = U[2, 1] * a + U[2, 2] * b
    end
    return psi
end

function two_gate!(psi, U, a, b)
    a == b && throw(ArgumentError("two-qubit gate needs distinct sites"))
    abit, bbit = 1 << (a - 1), 1 << (b - 1)
    for i in 0:(length(psi) - 1)
        (i & (abit | bbit)) == 0 || continue
        indices = (i + 1, (i | bbit) + 1, (i | abit) + 1, (i | abit | bbit) + 1)
        v = U * ComplexF64[psi[j] for j in indices]
        for j in 1:4
            psi[indices[j]] = v[j]
        end
    end
    return psi
end

function cz!(psi, a, b)
    bits = (1 << (a - 1)) | (1 << (b - 1))
    for i in 0:(length(psi) - 1)
        (i & bits) == bits && (psi[i + 1] = -psi[i + 1])
    end
    return psi
end

function t_gate!(psi, q)
    bit = 1 << (q - 1)
    phase = cis(pi / 4)
    for i in 0:(length(psi) - 1)
        (i & bit) != 0 && (psi[i + 1] *= phase)
    end
    return psi
end

function measure_z!(rng, psi, q)
    bit = 1 << (q - 1)
    p0 = sum(abs2(psi[i + 1]) for i in 0:(length(psi) - 1) if (i & bit) == 0)
    outcome = rand(rng) < p0 ? 0 : 1
    probability = outcome == 0 ? p0 : 1 - p0
    probability > 0 || error("sampled a zero-probability measurement outcome")
    for i in 0:(length(psi) - 1)
        ((i & bit) == 0) == (outcome == 0) || (psi[i + 1] = 0)
    end
    psi ./= sqrt(probability)
    return outcome
end

function measure!(rng, psi, q, basis)
    basis in (:X, :Z) || throw(ArgumentError("basis must be :X or :Z"))
    basis == :X && single_gate!(psi, H, q)
    outcome = measure_z!(rng, psi, q)
    basis == :X && single_gate!(psi, H, q)
    return outcome
end

function random_pair(rng, L)
    a = rand(rng, 1:L)
    b = rand(rng, 1:(L - 1))
    b >= a && (b += 1)
    return a, b
end

function paper_step!(rng, psi, L; p_cz=0.5, p_t=0.0, p_m=0.0, basis=:X)
    L >= 2 || throw(ArgumentError("L must be at least 2"))
    all(0 <= p <= 1 for p in (p_cz, p_t, p_m)) ||
        throw(ArgumentError("probabilities must lie in [0,1]"))
    rand(rng) < p_cz && cz!(psi, random_pair(rng, L)...)
    rand(rng) < p_t && t_gate!(psi, rand(rng, 1:L))
    if rand(rng) < p_m
        q = rand(rng, 1:L)
        outcome = measure!(rng, psi, q, basis)
        if basis == :Z
            outcome == 1 && single_gate!(psi, X, q)
            single_gate!(psi, H, q)
        end
        return (measured=true, site=q, outcome=outcome)
    end
    return (measured=false, site=0, outcome=0)
end

function ancilla_entropy(psi, L)
    length(psi) == 1 << (L + 1) || throw(ArgumentError("state has wrong dimension"))
    n = 1 << L
    p0 = sum(abs2, @view psi[1:n])
    p1 = sum(abs2, @view psi[(n + 1):(2n)])
    coherence = sum(psi[i] * conj(psi[n + i]) for i in 1:n)
    delta = sqrt(max(0.0, (p0 - p1)^2 + 4abs2(coherence)))
    eigenvalues = ((p0 + p1 + delta) / 2, (p0 + p1 - delta) / 2)
    return -sum(lambda * log(lambda) for lambda in eigenvalues if lambda > 1e-14)
end

function purification_trajectory(L; p_m, p_t=0.0, tf=2L^2,
                                 scramble_gates=round(Int, sqrt(10) * L),
                                 rng=MersenneTwister(1), cliffords=nothing)
    L >= 2 || throw(ArgumentError("L must be at least 2"))
    L <= 12 || throw(ArgumentError("dense reference is limited to L <= 12"))
    tf >= 0 || throw(ArgumentError("tf must be nonnegative"))
    cliffords === nothing && (cliffords = two_qubit_cliffords())
    psi = zeros(ComplexF64, 1 << (L + 1))
    psi[1:(1 << L)] .= inv(sqrt(1 << L))
    site = rand(rng, 1:L)
    two_gate!(psi, CNOT, site, L + 1)
    for _ in 1:scramble_gates
        a, b = random_pair(rng, L)
        two_gate!(psi, rand(rng, cliffords), a, b)
    end
    entropies = Float64[ancilla_entropy(psi, L)]
    for _ in 1:tf
        paper_step!(rng, psi, L; p_t, p_m, basis=:X)
        push!(entropies, ancilla_entropy(psi, L))
    end
    return entropies
end

function purification_ensemble(L; p_m, p_t=0.0, tf=2L^2, trajectories=16,
                               seed=1, scramble_gates=round(Int, sqrt(10) * L),
                               cliffords=nothing)
    trajectories > 0 || throw(ArgumentError("trajectories must be positive"))
    cliffords === nothing && (cliffords = two_qubit_cliffords())
    samples = Matrix{Float64}(undef, tf + 1, trajectories)
    for trial in 1:trajectories
        samples[:, trial] = purification_trajectory(L; p_m, p_t, tf, scramble_gates,
            rng=MersenneTwister(seed + trial - 1), cliffords)
    end
    means = vec(mean(samples; dims=2))
    sems = trajectories == 1 ? zeros(tf + 1) : vec(std(samples; dims=2)) ./ sqrt(trajectories)
    return (time=collect(0:tf), mean=means, sem=sems, samples=samples)
end

function fit_late_tail(time, mean_entropy, sem, L; min_points=4)
    length(time) == length(mean_entropy) == length(sem) ||
        throw(ArgumentError("time, mean, and sem must have equal length"))
    minimum_width = L^2 / 8
    best = nothing
    for i in eachindex(time), j in (i + min_points - 1):length(time)
        j <= length(time) || continue
        time[j] - time[i] >= minimum_width || continue
        mean_entropy[j] > 0 && sem[j] / mean_entropy[j] < 0.35 || continue
        ys = @view mean_entropy[i:j]
        all(y -> y > 0, ys) || continue
        xs = Float64.(time[i:j])
        logs = log.(ys)
        xbar, ybar = mean(xs), mean(logs)
        sxx = sum(abs2(x - xbar) for x in xs)
        sxx > 0 || continue
        slope = sum((xs[k] - xbar) * (logs[k] - ybar) for k in eachindex(xs)) / sxx
        slope < 0 || continue
        intercept = ybar - slope * xbar
        residual = sum(abs2(logs[k] - intercept - slope * xs[k]) for k in eachindex(xs))
        total = sum(abs2(y - ybar) for y in logs)
        r2 = total > 0 ? 1 - residual / total : 1.0
        candidate = (tau=-1 / slope, t_start=time[i], t_end=time[j], r2=r2)
        (best === nothing || candidate.r2 > best.r2) && (best = candidate)
    end
    return best
end

end
