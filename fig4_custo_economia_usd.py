#!/usr/bin/env python3
"""Regenerates Figure 4 (a and b) with monetary values in US dollars.

The original panels, extracted from the Word manuscript, were labelled in
Brazilian reais and in Portuguese. This script redraws them in USD (USD 1.00 =
BRL 5.11) using the same visual standard as the remaining figures of the paper
(serif type, navy/steel-blue palette, light grid, English labels).

Outputs (written to ./figs/):
    fig4a-ft-custo-operacional.png
    fig4b-economia-potencial.png
"""

from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent / "figs"
OUT.mkdir(exist_ok=True)

# House palette, sampled from the other figures of the manuscript.
NAVY = "#1a3a5c"
BLUE = "#2b6cb0"
FILL = "#c7e7f9"
GRID = "#e1e1e1"
TEXT = "#333333"

mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 13,
    "axes.edgecolor": TEXT,
    "axes.labelcolor": TEXT,
    "axes.linewidth": 0.9,
    "text.color": TEXT,
    "xtick.color": TEXT,
    "ytick.color": TEXT,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "legend.frameon": True,
    "legend.edgecolor": "#bbbbbb",
    "legend.framealpha": 1.0,
    "savefig.dpi": 412,
    "figure.dpi": 412,
})

# Weibull parameters (MLE) and cost of a corrective intervention.
BETA, ETA = 1.4449, 21.533
CF_USD = 1309.78          # BRL 6,692.97 at USD 1.00 = BRL 5.11


def fig4a():
    t = np.linspace(0.5, 60, 600)
    Ft = 1.0 - np.exp(-((t / ETA) ** BETA))
    cost = (BETA / ETA) * (t / ETA) ** (BETA - 1.0) * CF_USD

    fig, ax = plt.subplots(figsize=(10, 6.36))
    ax.set_axisbelow(True)
    ax.grid(True)

    l1, = ax.plot(t, Ft, color=NAVY, lw=2.6,
                  label=r"$F(t)$ — cumulative failure probability")
    ax.set_xlabel("Operating time without preventive intervention (days)")
    ax.set_ylabel(r"$F(t)$ — cumulative failure probability", color=NAVY)
    ax.tick_params(axis="y", colors=NAVY)
    ax.set_ylim(0, 1.02)
    ax.set_xlim(0, 61)

    ax2 = ax.twinx()
    l2, = ax2.plot(t, cost, color=BLUE, lw=2.6, ls="--",
                   label=r"Expected operating cost $h(t)\cdot C_f$")
    ax2.set_ylabel("Expected operating cost (USD/day)", color=BLUE)
    ax2.tick_params(axis="y", colors=BLUE)
    ax2.grid(False)

    ax.set_title(r"$F(t)$ and expected operating cost $h(t)\cdot C_f$", pad=12)
    ax.legend(handles=[l1, l2], loc="upper left", fontsize=12)

    fig.tight_layout()
    fig.savefig(OUT / "fig4a-ft-custo-operacional.png", facecolor="white")
    fig.savefig(OUT / "fig4a-ft-custo-operacional.pdf", facecolor="white")
    plt.close(fig)


def fig4b():
    # Annual cost per device (USD), Table V.
    rtf = 24483.34
    strategies = [("PMBC", 22044.37), ("RCM", 20124.85), ("Predictive", 14964.32)]

    x = np.arange(len(strategies))
    w = 0.36

    fig, ax = plt.subplots(figsize=(10, 6.36))
    ax.set_axisbelow(True)
    ax.grid(True, axis="y")

    ax.bar(x - w / 2, [rtf] * len(strategies), w, color=NAVY,
           edgecolor=NAVY, label="Current scenario (corrective/RTF)")
    ax.bar(x + w / 2, [c for _, c in strategies], w, color=FILL,
           edgecolor=BLUE, linewidth=1.2, label="Proposed scenario")

    for i, (_, cost) in enumerate(strategies):
        ax.text(i - w / 2, rtf + 350, f"USD {rtf:,.0f}",
                ha="center", va="bottom", fontsize=11.5)
        ax.text(i + w / 2, cost - 900, f"USD {cost:,.0f}",
                ha="center", va="top", fontsize=11.5, color=NAVY)
        saving = rtf - cost
        ax.text(i + w / 2, (cost + rtf) / 2,
                f"$-$USD {saving:,.0f}\n($-${saving / rtf * 100:.1f}%)",
                ha="center", va="center", fontsize=12,
                color=BLUE, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                          edgecolor=BLUE, linewidth=0.9))

    ax.set_xticks(x, [s for s, _ in strategies])
    ax.set_ylabel("Annual cost (USD/device/year)")
    ax.yaxis.set_major_formatter(
        mpl.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_ylim(0, rtf * 1.22)
    ax.set_title("Potential savings: current vs. proposed scenario", pad=12)
    ax.legend(loc="upper right", fontsize=12)

    fig.tight_layout()
    fig.savefig(OUT / "fig4b-economia-potencial.png", facecolor="white")
    fig.savefig(OUT / "fig4b-economia-potencial.pdf", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    fig4a()
    fig4b()
