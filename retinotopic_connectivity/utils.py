import numpy as np
from matplotlib import pyplot as plt

import numpy as np

def nature_style_plot(
    ax,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
    xticks=None,
    yticks=None,
    n_xticks=None,
    n_yticks=3,
    spine_width=2,
    tick_length=6,
    tick_width=2,
    fontsize=16,
    x_decimals=0,
    y_decimals=3,
    add_origin_padding=True,
    pad_fraction=0.02,
    format_xticklabels=True,
    format_yticklabels=True
):
    """
    Apply Nature-style formatting to matplotlib axes.

    Priority:
    1. explicit xticks / yticks
    2. n_xticks / n_yticks
    3. keep existing ticks
    """

    # Remove top/right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Style spines
    ax.spines["left"].set_linewidth(spine_width)
    ax.spines["bottom"].set_linewidth(spine_width)

    # Tick style
    ax.tick_params(
        axis="both",
        direction="out",
        length=tick_length,
        width=tick_width,
        labelsize=fontsize
    )

    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    # Set limits
    if xmin is not None and xmax is not None:
        ax.set_xlim(xmin, xmax)
    if ymin is not None and ymax is not None:
        ax.set_ylim(ymin, ymax)

    # Optional padding
    if add_origin_padding:
        xmin_current, xmax_current = ax.get_xlim()
        ymin_current, ymax_current = ax.get_ylim()

        xpad = pad_fraction * (xmax_current - xmin_current)
        ypad = pad_fraction * (ymax_current - ymin_current)

        ax.set_xlim(xmin_current - xpad, xmax_current)
        ax.set_ylim(ymin_current - ypad, ymax_current)

    # ----- X ticks -----
    if xticks is not None:
        ax.set_xticks(xticks)
        if format_xticklabels:
            ax.set_xticklabels([f"{t:.{x_decimals}f}" for t in xticks])
    elif n_xticks is not None:
        xmin_current, xmax_current = ax.get_xlim()
        if n_xticks == 3:
            ticks = [xmin_current, (xmin_current + xmax_current) / 2, xmax_current]
        elif n_xticks == 2:
            ticks = [xmin_current, xmax_current]
        else:
            ticks = np.linspace(xmin_current, xmax_current, n_xticks)

        ax.set_xticks(ticks)
        if format_xticklabels:
            ax.set_xticklabels([f"{t:.{x_decimals}f}" for t in ticks])

    # ----- Y ticks -----
    if yticks is not None:
        ax.set_yticks(yticks)
        if format_yticklabels:
            ax.set_yticklabels([f"{t:.{y_decimals}f}" for t in yticks])
    elif n_yticks is not None:
        ymin_current, ymax_current = ax.get_ylim()
        if n_yticks == 3:
            ticks = [ymin_current, (ymin_current + ymax_current) / 2, ymax_current]
        elif n_yticks == 2:
            ticks = [ymin_current, ymax_current]
        else:
            ticks = np.linspace(ymin_current, ymax_current, n_yticks)

        ax.set_yticks(ticks)
        if format_yticklabels:
            ax.set_yticklabels([f"{t:.{y_decimals}f}" for t in ticks])

    # Trim spines to final ticks
    xticks_final = ax.get_xticks()
    yticks_final = ax.get_yticks()

    if len(xticks_final) > 0:
        ax.spines["bottom"].set_bounds(xticks_final[0], xticks_final[-1])

    if len(yticks_final) > 0:
        ax.spines["left"].set_bounds(yticks_final[0], yticks_final[-1])

    return ax

