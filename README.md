# Low-Rank Stabilizer State Simulator

[![Tests](https://github.com/ka784-jpg/Low-Rank-Stabilizer-State-Simulator/actions/workflows/julia-tests.yml/badge.svg)](https://github.com/ka784-jpg/Low-Rank-Stabilizer-State-Simulator/actions/workflows/julia-tests.yml)
[![Code license: MIT](https://img.shields.io/badge/code_license-MIT-blue.svg)](LICENSE)

Julia simulation kernels and example drivers for low-magic quantum dynamics.
The repository includes low-rank stabilizer-decomposition (LRSD) evolution,
Bell-sampling stabilizer-nullity calculations, and a small-system dense
reference for ancilla purification, with tests and published figure previews.

**Paper:** Kemal Aziz, Haining Pan, Michael J. Gullans, and J. H. Pixley,
[*Classical Simulations of Low Magic Quantum Dynamics*](https://arxiv.org/abs/2508.20252).
Figure references below follow [arXiv version 3](https://arxiv.org/abs/2508.20252v3).

[Run an example](#quickstart) | [Methods](#methods) |
[Paper figures and drivers](scripts/README.md) |
[Reproducibility](paper/REPRODUCIBILITY.md) | [Figure-to-code map](paper/FIGURE_MAP.md)

![Published Figure 4: LRSD representation size versus system size, truncation cutoff, and T-gate rate in the X-measurement model.](paper/figures/fig04.png)

**What the calculation probes.** How large a classical representation is
needed as a monitored circuit accumulates non-Clifford resources? Figure 4
compares the LRSD representation size across system sizes, truncation
thresholds, and T-gate rates. It shows the paper's results, not output from
the small example below. [Panel context and simulation driver](scripts/README.md#lrsd-truncation-and-cost)
and [figure attribution](paper/figures/README.md).

## Explore the results

| LRSD truncation and cost | Ancilla purification | Bell sampling and magic |
| --- | --- | --- |
| [![Published Figure 3: truncation diagnostics.](paper/figures/fig03.png)](scripts/README.md#lrsd-truncation-and-cost) | [![Published Figure 5: purification transition.](paper/figures/fig05.png)](scripts/README.md#ancilla-purification) | [![Published Figure 6: magic scaling.](paper/figures/fig06.png)](scripts/README.md#bell-sampling-and-magic) |
| [Figures 3, 4, 17 and driver](scripts/README.md#lrsd-truncation-and-cost) | [Figures 5, 14-16 and driver](scripts/README.md#ancilla-purification) | [Figures 6, 8 and driver](scripts/README.md#bell-sampling-and-magic) |

The [complete figure map](paper/FIGURE_MAP.md) identifies which calculations
are executable and which still lack an implementation or production data.
The dense ancilla reference is limited to small systems; it is not the
scalable LRSD engine used for the paper's large-system purification results.

## Quickstart

Install [Julia](https://julialang.org/downloads/) (tested with Julia 1.12),
then run these commands from a terminal:

```sh
git clone https://github.com/ka784-jpg/Low-Rank-Stabilizer-State-Simulator.git
cd Low-Rank-Stabilizer-State-Simulator
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/fig34_lrsd.jl --L 4,6 --pm 0.6,0.8 --basis X --eta 1 --beta 1 --cutoff 1e-15 --trajectories 4 --seed 2026 --out outputs/quickstart.csv
```

This small LRSD run writes 16 trajectory rows with the parameters, seed,
logical-term count, representation entries, and runtime. It samples the
observable in Figure 4(a), but does not recreate that figure's large-system
statistics. The first Julia run includes compilation.

To plot the CSV with Python and Matplotlib:

```sh
python -m pip install matplotlib
python scripts/plot_examples.py fig4a outputs/quickstart.csv outputs/quickstart.png
```

Use a fresh output filename to retain previous results: the existing drivers
can overwrite their CSV destination. No cluster, scheduler, or credentials are
needed. [More examples](#small-executable-examples) cover truncation sweeps,
purification, and nullity.

## Methods

- `LRSD Scaling Analysis`: original low-rank stabilizer-decomposition (LRSD) kernel, with a seedable CSV driver for logical-term counts and representation size. The driver covers Fig. 3(a), Fig. 4, and trajectory-level input for Fig. 17.
- `Stabilizer Nullity (Magic) Computation`: original Z-basis Clifford-plus-terminal-T and Bell-sampling nullity kernel. The driver covers Fig. 6(b,c) data, but not the finite-size transition-location analysis of Fig. 6(a).
- `src/PaperAllToAll.jl`: new **dense, small-system reference** for the paper's ancilla purification protocol. It evolves the system plus one reference qubit exactly, with all-to-all CZ/T/X-measurement steps, a random two-qubit Clifford scrambling stage, ancilla entropy, and late-time exponential-window fitting. It supports Fig. 5 and Figs. 14-16 at small L; it is not the scalable LRSD ancilla implementation used for the published L=32-128 results.

The Z-basis circuit measures Z and resets the measured qubit to `|+>` as specified in the paper. T and CZ commute, which permits the LRSD driver to apply T before its CZ/measurement layer without changing the per-step circuit.

## Install and test

Install [Julia](https://julialang.org/downloads/) (tested locally with Julia 1.12) and run from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
```

The first run compiles the original kernels and can take a minute. The optional plot helper requires Python and Matplotlib (`python -m pip install matplotlib`). No HPC scheduler or campaign is invoked by these examples.

## Small executable examples

```sh
# Figure 3(a): Z-basis term count vs cutoff; common random seed per cutoff.
julia --project=. scripts/fig34_lrsd.jl --L 4 --pm 0.6,0.8 --basis Z --eta 1 --beta 1 --cutoff 0.0001,0.01,0.1 --trajectories 8 --out outputs/fig3a.csv

# Figure 4(a): X-basis LRSD representation size vs L.
julia --project=. scripts/fig34_lrsd.jl --L 4,6,8 --pm 0.6,0.7,0.8,0.9 --basis X --eta 1 --beta 1 --cutoff 1e-15 --trajectories 8 --out outputs/fig4.csv

# Figure 5 / Appendix Fig. 14: exact small-L ancilla purification.
julia --project=. scripts/fig5_purification.jl --L 4,6 --pm 0.2,0.26,0.32 --eta 0 --trajectories 16 --out outputs/fig5.csv

# Figure 6(b,c): Z-basis Bell-sampling stabilizer nullity.
julia --project=. scripts/fig6_bell_magic.jl --L 16,24 --pm 0.5 --eta 12 --beta 0.9,1.0 --trajectories 8 --out outputs/fig6.csv
```

The `--beta` argument accepts a comma-separated list in the Fig. 6 driver. For every point, `pT = eta/L^beta`; requests with `pT > 1` are rejected. Increase sizes and trajectories deliberately: LRSD can grow exponentially with magic, and the dense reference uses `2^(L+1)` amplitudes and is limited to `L <= 12`.

```sh
python scripts/plot_examples.py fig3a outputs/fig3a.csv outputs/fig3a.png
python scripts/plot_examples.py fig4 outputs/fig4.csv outputs/fig4.png
python scripts/plot_examples.py fig14 outputs/fig5.csv outputs/fig14.png
python scripts/plot_examples.py fig6b outputs/fig6.csv outputs/fig6b.png
python scripts/plot_examples.py fig6c outputs/fig6.csv outputs/fig6c.png
```

For Fig. 4(b), vary `--cutoff` at fixed `--pm` and `--eta` and plot with `fig4b`. For Fig. 4(c), vary `--eta` at fixed `--pm` and `--cutoff` and plot with `fig4c`. The Fig. 17 trajectory histogram uses `fig17` on X-basis LRSD output. The Fig. 6(c) helper plots `M/log(L)` from the same samples as Fig. 6(b).

CSV output contains the parameter values and random seed for each parameter point. These are **new simulations**; published raw data, uncertainty analyses, phase-boundary fits, and high-statistics figures are not bundled. The figure map also records important algorithms not yet available here.

## License

Repository-authored code and documentation are released under the [MIT License](LICENSE). The [published figure previews](paper/figures/README.md) retain the paper's CC BY 4.0 license, not MIT. Julia and Python dependencies are obtained separately and retain their own licenses; no third-party source is vendored in this repository.
