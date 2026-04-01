import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import FormatStrFormatter

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import FormatStrFormatter

def plot_dva_bar(shell_vals):
    """
    Bar plot of discrete DVA-distance connectivity with SEM error bars.

    shell_vals : dict
        shell_vals[k] = array of all M[i,j] such that |i-j| = k
    """

    labels = [
        "same DVA", "±1 DVA", "±2 DVA", "±3 DVA",
        "±4 DVA", "±5 DVA", "≥6 DVA"
    ]

    # Means and SEMs for first six shells
    means = [shell_vals[k].mean() for k in range(6)]
    sems  = [shell_vals[k].std() / np.sqrt(len(shell_vals[k])) for k in range(6)]

    # Pool distant bins (k >= 6)
    far_vals = np.concatenate([v for k, v in shell_vals.items() if k >= 6])
    means.append(far_vals.mean())
    sems.append(far_vals.std() / np.sqrt(len(far_vals)))

    # Truncated hot colormap (avoid white)
    base_cmap = cm.hot
    trunc_cmap = colors.LinearSegmentedColormap.from_list(
        "hot_truncated",
        base_cmap(np.linspace(0.0, 0.85, 256))
    )

    norm = colors.Normalize(vmin=min(means), vmax=max(means))
    bar_colors = list(trunc_cmap(norm(means)))

    # Lighten first two bars
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

    # Plot
    plt.figure()
    plt.bar(
        range(len(means)),
        means,
        yerr=sems,
        color=bar_colors,
        capsize=0,
        ecolor="0.5",
        error_kw=dict(linewidth=1)
    )

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ymin = 0
    ymax = max(np.array(means) + np.array(sems))
    ax.set_ylim(ymin, ymax * 1.05)
    ax.set_yticks([ymin, ymax / 2, ymax])
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    plt.xticks(range(len(means)), labels, rotation=30)
    plt.ylabel("Mean streamline density")

    plt.tight_layout()
    plt.show()


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def plot_radial_shells(M, out_png=None):
    """
    Plot mean connectivity as a function of distance from the diagonal |i-j|.
    """
    n = M.shape[0]

    ks = []
    means = []

    for k in range(n):
        vals = []
        for i in range(n):
            for j in range(n):
                if abs(i - j) == k:
                    vals.append(M[i, j])
        if len(vals) > 0:
            ks.append(k)
            means.append(np.mean(vals))

    ks = np.array(ks)
    means = np.array(means)

    plt.figure()
    plt.plot(ks, means, marker="o")

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ymin = 0
    ymax = means.max()
    ax.set_ylim(ymin, ymax * 1.05)
    ax.set_yticks([ymin, ymax / 2, ymax])
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    plt.xlabel("DVA distance (|ΔDVA|)")
    plt.ylabel("Mean streamline density")

    plt.tight_layout()

    if out_png is not None:
        plt.savefig(out_png, dpi=300)
        plt.close()
    else:
        plt.show()


import numpy as np

def compute_shell_vals(M):
    """
    Compute shell values from a connectivity matrix.

    Parameters
    ----------
    M : (N, N) array
        Connectivity matrix.

    Returns
    -------
    shell_vals : dict
        shell_vals[k] = 1D array of all M[i,j] with |i-j| = k
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

def shell_stats(shell_vals):
    """
    Compute mean and SEM for each shell.

    Returns
    -------
    means : dict
    sems  : dict
    """
    means = {}
    sems  = {}

    for k, vals in shell_vals.items():
        means[k] = vals.mean()
        sems[k]  = vals.std() / np.sqrt(len(vals))

    return means, sems

def dva_groups(shell_vals, max_explicit=6):
    """
    Group shells into discrete DVA bins:
    0, 1, 2, ..., max_explicit-1, and >= max_explicit
    """
    means = []
    sems  = []

    # explicit shells
    for k in range(max_explicit):
        vals = shell_vals[k]
        means.append(vals.mean())
        sems.append(vals.std() / np.sqrt(len(vals)))

    # pooled distant shells
    far_vals = np.concatenate(
        [v for k, v in shell_vals.items() if k >= max_explicit]
    )
    means.append(far_vals.mean())
    sems.append(far_vals.std() / np.sqrt(len(far_vals)))

    return np.array(means), np.array(sems)


def mat_to_barplot(M):
    """
    Wrapper: from connectivity matrix to DVA bar plot with SEM.
    """
    shell_vals = compute_shell_vals(M)
    plot_dva_bar(shell_vals)


def mat_to_radial_shells(M):
    """
    Wrapper: from connectivity matrix to radial shell plot.
    """
    plot_radial_shells(M)

def mat_to_connectivity_summary(M):
    """
    Generate both DVA bar plot and radial shell curve from a matrix.
    """
    shell_vals = compute_shell_vals(M)
    plot_dva_bar(shell_vals)
    plot_radial_shells(M)
