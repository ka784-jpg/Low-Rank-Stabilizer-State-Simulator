# Reproducibility and scope

[Repository overview](../README.md) | [Run an example](../README.md#quickstart) |
[Visual results guide](../scripts/README.md) | [Figure-to-code map](FIGURE_MAP.md)

## What you can run

- The original LRSD kernel, through `scripts/fig34_lrsd.jl`, records logical-term
  counts and representation size for small all-to-all monitored circuits.
- `scripts/fig6_bell_magic.jl` samples stabilizer nullity in the paper's Z-basis
  terminal-T construction. This is not the cutoff-dependent nullity analysis
  of Figure 3(b).
- `scripts/fig5_purification.jl` uses a dense reference for system-plus-ancilla
  trajectories and exponential-window fitting. It is limited to `L <= 12` and
  is not the large-system LRSD purification implementation.

The [main-page commands](../README.md#small-executable-examples) give small
parameter sets; the [visual guide](../scripts/README.md) shows the corresponding
published figures. The root `Project.toml` and `Manifest.toml` describe the Julia
environment. Run `julia --project=. test/runtests.jl` for the supplied tests.

## What the published images mean

The images are previews of Figures 1-19 from
[*Classical Simulations of Low Magic Quantum Dynamics*, arXiv:2508.20252v3](https://arxiv.org/abs/2508.20252v3),
by Kemal Aziz, Haining Pan, Michael J. Gullans, and J. H. Pixley.
They were rasterized from the figure PDFs in that version's source archive.
The [asset list](figures/README.md) records their source filenames and CC BY 4.0
terms. They are not outputs of the bundled small examples.

The code is MIT-licensed; that does not replace the figure license or the
licenses of separately installed dependencies.

## What is not reproduced here

The repository does not bundle the paper's large raw datasets, production-scale
ensembles, complete phase-boundary/exponent fitting, or all appendix algorithms.
The [coverage table](FIGURE_MAP.md) records these gaps figure by figure.
Its absent scalable ancilla workflow, cutoff-dependent nullity driver, and
several graph/entanglement analyses should not be inferred from the presence
of a published image.

A fixed seed supports rerunning a point in the same implementation, not a
cross-version or cross-platform bitwise guarantee. Resetting the seed across
LRSD cutoffs does not guarantee the same circuit record if control flow changes
random-number consumption. Exact pairing would require explicit record replay.

CSV files retain per-point seeds and parameters. Runtime columns are inherently
machine-dependent. The current drivers can overwrite the requested output
file, so use distinct output paths when comparing runs. Tests and executable
examples check implementation behavior; they do not by themselves validate
the published critical exponents or large-system claims.
