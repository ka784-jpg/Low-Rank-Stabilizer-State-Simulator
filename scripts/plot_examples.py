"""Plot CSVs produced by the paper example drivers (not the published data)."""

import argparse
import csv
import math
from collections import defaultdict
from statistics import mean

import matplotlib.pyplot as plt


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def aggregate(rows, keys, value):
    grouped = defaultdict(list)
    for row in rows:
        grouped[tuple(row[key] for key in keys)].append(float(row[value]))
    return {key: mean(values) for key, values in grouped.items()}


def plot_fig3a(ax, rows):
    if len({(r["L"], r["eta"], r["beta"], r["basis"]) for r in rows}) != 1:
        raise ValueError("Fig. 3(a) view needs one L/eta/beta/basis")
    if rows[0]["basis"] != "Z":
        raise ValueError("Fig. 3(a) is the Z-basis model")
    values = aggregate(rows, ("pm", "cutoff"), "logical_terms")
    for pm in sorted({key[0] for key in values}, key=float):
        points = sorted(
            ((float(cutoff), result) for (rate, cutoff), result in values.items()
             if rate == pm and float(cutoff) > 0),
            key=lambda point: point[0],
        )
        if points:
            ax.plot(*zip(*points), marker="o", label=f"$p_m={pm}$")
    ax.set(xscale="log", xlabel=r"Cutoff $\epsilon$", ylabel=r"Mean logical terms $K$")
    ax.legend()


def plot_fig4a(ax, rows):
    combos = {(r["cutoff"], r["eta"], r["beta"], r["basis"]) for r in rows}
    if len(combos) != 1 or next(iter(combos))[3] != "X":
        raise ValueError("Fig. 4 view needs a single cutoff/eta/beta and X basis")
    values = aggregate(rows, ("pm", "L"), "representation_entries")
    for pm in sorted({key[0] for key in values}, key=float):
        points = sorted(((int(L), result) for (rate, L), result in values.items()
                         if rate == pm), key=lambda point: point[0])
        ax.plot(*zip(*points), marker="o", label=f"$p_m={pm}$")
    ax.set(xscale="log", yscale="log", xlabel=r"System size $L$",
           ylabel=r"Mean representation entries $\mathcal{N}$")
    ax.legend()


def plot_fig4b(ax, rows):
    if len({(r["pm"], r["eta"], r["beta"], r["basis"]) for r in rows}) != 1:
        raise ValueError("Fig. 4(b) view needs one pm/eta/beta/basis")
    if rows[0]["basis"] != "X":
        raise ValueError("Fig. 4(b) is the X-basis model")
    values = aggregate(rows, ("cutoff", "L"), "representation_entries")
    for cutoff in sorted({key[0] for key in values}, key=float):
        points = sorted(((int(L), result) for (value, L), result in values.items()
                         if value == cutoff), key=lambda point: point[0])
        ax.plot(*zip(*points), marker="o", label=rf"$\epsilon={cutoff}$")
    ax.set(xscale="log", yscale="log", xlabel=r"System size $L$",
           ylabel=r"Mean representation entries $\mathcal{N}$")
    ax.legend()


def plot_fig4c(ax, rows):
    if len({(r["pm"], r["cutoff"], r["beta"], r["basis"]) for r in rows}) != 1:
        raise ValueError("Fig. 4(c) view needs one pm/cutoff/beta/basis")
    if rows[0]["basis"] != "X":
        raise ValueError("Fig. 4(c) is the X-basis model")
    values = aggregate(rows, ("eta", "L"), "representation_entries")
    for eta in sorted({key[0] for key in values}, key=float):
        points = sorted(((int(L), result) for (value, L), result in values.items()
                         if value == eta), key=lambda point: point[0])
        ax.plot(*zip(*points), marker="o", label=rf"$\eta={eta}$")
    ax.set(xscale="log", yscale="log", xlabel=r"System size $L$",
           ylabel=r"Mean representation entries $\mathcal{N}$")
    ax.legend()


def plot_fig5(ax, rows):
    points = {}
    for row in rows:
        if row["tau"]:
            points[(row["pm"], int(row["L"]))] = float(row["tau"])
    if not points:
        raise ValueError("No late-tail fits passed the paper's window criteria")
    for pm in sorted({key[0] for key in points}, key=float):
        series = sorted(((L, tau) for (rate, L), tau in points.items() if rate == pm))
        ax.plot(*zip(*series), marker="o", label=f"$p_m={pm}$")
    ax.set(xlabel=r"System size $L$", ylabel=r"Decay time $\tau$")
    ax.legend()


def plot_fig14(ax, rows):
    combos = {(r["L"], r["eta"], r["beta"]) for r in rows}
    if len(combos) != 1:
        raise ValueError("Time-series view needs a single L/eta/beta")
    for pm in sorted({r["pm"] for r in rows}, key=float):
        series = sorted(((int(r["t"]), float(r["mean_SQ"])) for r in rows
                         if r["pm"] == pm), key=lambda point: point[0])
        ax.plot(*zip(*series), label=f"$p_m={pm}$")
    ax.set(xlabel=r"Time $t$", ylabel=r"Mean ancilla entropy $S_Q$")
    ax.legend()


def plot_fig6b(ax, rows):
    pms = {r["pm"] for r in rows}
    if len(pms) != 1:
        raise ValueError("Fig. 6(b) view needs a single measurement probability")
    values = aggregate(rows, ("beta", "L"), "M_estimate")
    for beta in sorted({key[0] for key in values}, key=float):
        series = sorted(((int(L), result) for (rate, L), result in values.items()
                         if rate == beta), key=lambda point: point[0])
        ax.plot(*zip(*series), marker="o", label=rf"$\beta={beta}$")
    ax.set(xlabel=r"System size $L$", ylabel=r"Mean sampled nullity $M$")
    ax.legend()


def plot_fig6c(ax, rows):
    pms = {r["pm"] for r in rows}
    if len(pms) != 1:
        raise ValueError("Fig. 6(c) view needs one measurement probability")
    values = aggregate(rows, ("L", "beta"), "M_estimate")
    for L in sorted({key[0] for key in values}, key=int):
        series = sorted(((float(beta), result / math.log(int(L)))
                         for (size, beta), result in values.items() if size == L),
                        key=lambda point: point[0])
        ax.plot(*zip(*series), marker="o", label=rf"$L={L}$")
    ax.set(xlabel=r"$\beta$", ylabel=r"Mean sampled nullity $M/\log L$")
    ax.legend()


def plot_fig17(ax, rows):
    if len({(r["L"], r["pm"], r["beta"], r["cutoff"], r["basis"])
            for r in rows}) != 1:
        raise ValueError("Fig. 17 view needs one L/pm/beta/cutoff/basis")
    if rows[0]["basis"] != "X":
        raise ValueError("Fig. 17 is the X-basis model")
    for eta in sorted({r["eta"] for r in rows}, key=float):
        samples = [float(r["logical_terms"]) for r in rows if r["eta"] == eta]
        ax.hist(samples, bins=min(30, max(5, len(set(samples)))), density=True,
                histtype="step", linewidth=1.5, label=rf"$\eta={eta}$")
    ax.set(xlabel=r"Logical terms $K$", ylabel="Density")
    ax.legend()


PLOTTERS = {"fig3a": plot_fig3a, "fig4": plot_fig4a, "fig4a": plot_fig4a,
            "fig4b": plot_fig4b, "fig4c": plot_fig4c, "fig5": plot_fig5,
            "fig14": plot_fig14, "fig6b": plot_fig6b, "fig6c": plot_fig6c,
            "fig17": plot_fig17}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("figure", choices=sorted(PLOTTERS))
    parser.add_argument("csv")
    parser.add_argument("output")
    args = parser.parse_args()
    rows = read_rows(args.csv)
    if not rows:
        raise ValueError("Input CSV is empty")
    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    PLOTTERS[args.figure](ax, rows)
    ax.grid(alpha=0.2)
    fig.savefig(args.output, dpi=200)
    print(args.output)


if __name__ == "__main__":
    main()
