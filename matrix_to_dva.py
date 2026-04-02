#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import FormatStrFormatter
from utils import nature_style_plot

# ============================================================
# Utilities
# ============================================================

def compute_shell_vals(M):
    """
    shell_vals[k] = array of all M[i,j] with |i-j| = k
    """
    n = M.shape[0]
    shell_vals = {}

    for k in range(n):
        vals = []
        for i in range(n):
            for j in range(n):
                if abs(i - j) == k:
                    vals.append(M[i, j])
        shell_vals[k] = np.array(vals)

    return shell_vals


# ============================================================
# Bar plot
# ============================================================

def plot_dva_bar(
    shell_vals,
    out_png: Path,
    bar_x_dim: float = 85,
    base_width: float = 6,
    base_height: float = 4,
    y_lim=None,
    fontsize=12,
    y_decimals=3
    ):
    """
    Bar plot of discrete DVA-distance connectivity with SEM error bars.
    """

    labels = [
        "same DVA", "±1 DVA", "±2 DVA", "±3 DVA",
        "±4 DVA", "±5 DVA", "≥6 DVA"
    ]

    # --------------------------------------------------------
    # Compute means + SEMs
    # --------------------------------------------------------
    means = [shell_vals[k].mean() for k in range(6)]
    sems  = [shell_vals[k].std() / np.sqrt(len(shell_vals[k])) for k in range(6)]

    far_vals = np.concatenate([v for k, v in shell_vals.items() if k >= 6])
    means.append(far_vals.mean())
    sems.append(far_vals.std() / np.sqrt(len(far_vals)))

    # --------------------------------------------------------
    # Shrink control
    # --------------------------------------------------------
    shrink = bar_x_dim / 100.0
    shrink = np.clip(shrink, 0.6, 1.2)

    fig_w = base_width * (0.85 + 0.15 * shrink)
    fig_h = base_height
    bar_w = 0.65
    x = np.arange(len(means)) * shrink

    # --------------------------------------------------------
    # Colormap
    # --------------------------------------------------------
    base_cmap = cm.hot
    trunc_cmap = colors.LinearSegmentedColormap.from_list(
        "hot_truncated",
        base_cmap(np.linspace(0.0, 0.85, 256))
    )

    norm = colors.Normalize(vmin=min(means), vmax=max(means))
    bar_colors = list(trunc_cmap(norm(means)))

    def lighten(color, amount=0.35):
        r, g, b, a = color
        return (
            r + (1 - r) * amount,
            g + (1 - g) * amount,
            b + (1 - b) * amount,
            a,
        )

    bar_colors[0] = lighten(bar_colors[0])
    bar_colors[1] = lighten(bar_colors[1])

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------
    plt.figure(figsize=(fig_w, fig_h))

    plt.bar(
        x,
        means,
        width=bar_w,
        yerr=sems,
        color=bar_colors,
        capsize=0,
        ecolor="0.5",
        error_kw=dict(linewidth=1),
    )

    ax = plt.gca()

    if y_lim is not None:
        ymin, ymax = y_lim
    else:
        ymin = 0
        ymax = max(np.array(means) + np.array(sems)) * 1.05

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30)

    nature_style_plot(
        ax,
        ymin=ymin,
        ymax=ymax,
        xticks=x,
        yticks=[ymin, (ymin + ymax) / 2, ymax],
        fontsize=fontsize,
        y_decimals=y_decimals,
        format_xticklabels=False,
        add_origin_padding=False,
    )

    plt.ylabel("Mean streamline density")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"✓ Saved bar plot → {out_png}")


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate DVA bar plot from a connectivity matrix CSV"
    )
    parser.add_argument(
        "matrix_csv",
        type=Path,
        help="CSV file containing the connectivity matrix"
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output PNG file (default: <matrix>_dva_barplot.png)"
    )
    parser.add_argument(
        "--bar-x-dim",
        type=float,
        default=85,
        help="Horizontal compression factor (default: 85)"
    )
    parser.add_argument(
        "--base-width",
        type=float,
        default=6,
        help="Base figure width (default: 6)"
    )
    parser.add_argument(
        "--base-height",
        type=float,
        default=4,
        help="Base figure height (default: 4)"
    )
    parser.add_argument(
        "--y-lim",
        type=str,
        default=None,
        help="Y-axis limits as min,max (example: 0,0.030)"
    )
    parser.add_argument(
        "--y-decim",
        type=str,
        default=None,
        help="Y-axis decimals (defaul 3)"
    )

    args = parser.parse_args()

    # Parse y-limits
    y_lim = None
    if args.y_lim is not None:
        try:
            y_min, y_max = map(float, args.y_lim.split(","))
            y_lim = (y_min, y_max)
        except ValueError:
            raise ValueError("Invalid --y-lim format. Use: min,max (example: 0,0.030)")

    M = np.loadtxt(args.matrix_csv, delimiter=",")
    if args.y_decim is not None:
        y_decim=args.y_decim
    else:
        y_decim=3
    out_png = (
        args.out
        if args.out is not None
        else args.matrix_csv.with_suffix("").with_name(
            args.matrix_csv.stem + "_dva_barplot.png"
        )
    )

    shell_vals = compute_shell_vals(M)

    plot_dva_bar(
        shell_vals,
        out_png=out_png,
        bar_x_dim=args.bar_x_dim,
        base_width=args.base_width,
        base_height=args.base_height,
        y_lim=y_lim,
        y_decimals=y_decim
    )


if __name__ == "__main__":
    main()
