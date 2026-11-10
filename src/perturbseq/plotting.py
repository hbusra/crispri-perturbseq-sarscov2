"""Shared palette and axis style.

A three-slot categorical set plus a blue-red diverging pair with a grey
neutral, checked for colour-vision deficiency separation on an all-pairs
basis. Every series is also direct-labelled or legended.
"""
from __future__ import annotations

import matplotlib.pyplot as plt

# ----------------------------------------------------------------- palette
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RED = "#e34948"
NEUTRAL = "#b4b3ad"
GRID = "#e6e5e1"
INK = "#0b0b0b"
INK2 = "#52514e"

STATE_COLOR = {"uninfected": BLUE, "infected": ORANGE, "borderline": NEUTRAL}
STATE_LABEL = {"uninfected": "bystander", "infected": "infected", "borderline": "borderline"}

plt.rcParams.update({
    "figure.dpi": 160,
    "savefig.dpi": 160,
    "savefig.bbox": "tight",
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.5,
    "axes.titleweight": "semibold",
    "axes.edgecolor": "#c9c8c3",
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK2,
    "text.color": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.frameon": False,
    "legend.fontsize": 7.5,
    "lines.linewidth": 2.0,
})


def style(ax, grid_axis="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
