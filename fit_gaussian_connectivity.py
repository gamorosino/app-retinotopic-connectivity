#!/usr/bin/env python3
"""
Fit Gaussian models to each row of a connectivity matrix.

Can be run standalone:
    python fit_gaussian_connectivity.py --matrix M.csv --ecc-bins 0_1,1_2,2_3,...

Or imported:
    from fit_gaussian_connectivity import fit_matrix_rows, ecc_bin_centers
"""

import numpy as np
import argparse
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# === Utils
# ============================================================

def save_figure(path, dpi=300):
    """Save plot in PNG and vector formats."""
    plt.savefig(path.with_suffix(".png"), dpi=dpi)
    plt.savefig(path.with_suffix(".pdf"))
    plt.savefig(path.with_suffix(".svg"))


# ============================================================
# === Gaussian Model Definitions
# ============================================================

def gaussian(x, k, mu, sigma, c):
    """Full Gaussian with amplitude, mean, width, baseline."""
    return k * np.exp(-0.5 * ((x - mu) / sigma)**2) + c


def gaussian_fixed_mu(x, k, sigma, c, mu_fixed):
    """Gaussian centered at a fixed mu."""
    return k * np.exp(-0.5 * ((x - mu_fixed) / sigma)**2) + c

def truncated_gaussian(x, A, mu, sigma, c, xmin, xmax, normalize=False):
    """
    Gaussian evaluated only in [xmin, xmax].
    If normalize=True → rescale so that sum(model) == sum(y_obs).
    """
    g = A * np.exp(-0.5 * ((x - mu) / sigma)**2)

    # truncate outside
    mask = (x >= xmin) & (x <= xmax)
    g = g * mask

    # optional normalization
    if normalize:
        total = g.sum()
        if total > 0:
            g = g / total  # becomes a probability mass
    return g + c

# ============================================================
# === Utilities
# ============================================================

def ecc_bin_centers(ecc_bins):
    """Convert ecc bin strings like '1_2' to numeric centers."""
    centers = []
    for b in ecc_bins:
        lo, hi = map(float, b.split("_"))
        centers.append((lo + hi) / 2)
    return np.array(centers)


# ============================================================
# === Gaussian Fitting for a Single Row
# ============================================================

def fit_row_gaussian(y, x, mu_init, fix_mu=False):
    """
    Fit a single row y = D[i,:] with a Gaussian.
    x = numeric centers of ecc bins
    mu_init = diagonal center (x[i])
    """
    y = np.asarray(y, dtype=float)

    # Initial guesses
    k0 = max(y) - min(y)
    c0 = min(y)
    sigma0 = (x.max() - x.min()) / 4

    if fix_mu:
        # Fit only k, sigma, c
        def model(x_, k, sigma, c):
            return gaussian_fixed_mu(x_, k, sigma, c, mu_fixed=mu_init)

        p0 = [k0, sigma0, c0]
        bounds = (
            [0.0, 1e-3, -np.inf],
            [np.inf, np.inf, np.inf]
        )

    else:
        # Fit full k, mu, sigma, c
        def model(x_, k, mu, sigma, c):
            return gaussian(x_, k, mu, sigma, c)

        mu_min = x.min()
        mu_max = x.max()

        p0 = [k0, mu_init, sigma0, c0]
        bounds = (
            [0.0, mu_min, 1e-3, -np.inf],
            [np.inf, mu_max, np.inf, np.inf]
        )

    popt, pcov = curve_fit(model, x, y, p0=p0, bounds=bounds, maxfev=20000)
    yfit = model(x, *popt)
    return popt, yfit


def fit_row_truncated_gaussian(y, x, mu_init, fix_mu=False, normalize=False):
    """
    Fit truncated Gaussian to a row y.
    - y: 1×N data row
    - x: centers (or indices)
    - mu_init: diagonal location
    - fix_mu: True → μ = mu_init
    - normalize: True → normalize Gaussian after truncation
    """

    y = np.asarray(y, float)
    xmin, xmax = x.min(), x.max()

    # initial parameters
    A0 = y.max() - y.min()
    c0 = y.min()
    sigma0 = (x.max()-x.min()) / 4

    if fix_mu:
        def model(x_, A, sigma, c):
            return truncated_gaussian(x_, A, mu_init, sigma, c,
                                      xmin, xmax, normalize=normalize)

        p0 = [A0, sigma0, c0]
        bounds = ([0, 1e-3, -np.inf], [np.inf, np.inf, np.inf])

    else:
        def model(x_, A, mu, sigma, c):
            return truncated_gaussian(x_, A, mu, sigma, c,
                                      xmin, xmax, normalize=normalize)

        p0 = [A0, mu_init, sigma0, c0]
        bounds = (
            [0, xmin, 1e-3, -np.inf],
            [np.inf, xmax, np.inf, np.inf]
        )

    popt, _ = curve_fit(model, x, y, p0=p0, bounds=bounds, maxfev=30000)
    yfit = model(x, *popt)

    return popt, yfit


# ============================================================
# === Fit All Rows of Matrix
# ============================================================

def fit_matrix_rows(D, ecc_bins, fix_mu=False):
    """
    Fit every row of D with a Gaussian.
    Returns:
        params  (N × P) matrix of Gaussian parameters
        fitted  (N × N) matrix of reconstructed Gaussian rows
    """
    x = ecc_bin_centers(ecc_bins)
    n = D.shape[0]

    params = []
    fitted = np.zeros_like(D)

    for i in range(n):
        row = D[i, :]
        mu_init = x[i]  # diagonal
        popt, yfit = fit_row_gaussian(row, x, mu_init, fix_mu=fix_mu)
        params.append(popt)
        fitted[i, :] = yfit

    return np.array(params), fitted


def fit_matrix_rows_fixed_mu(D, ecc_bins):
    """
    Fit each row i of matrix D with a Gaussian centered at μ=i.
    Free parameters: A (amplitude), sigma (spread).

    Returns:
        params_fixed: array of shape (N, 3) with columns [A_i, sigma_i, mu_i]
        fitted_fixed: N×N fitted matrix
    """

    from scipy.optimize import curve_fit
    import numpy as np

    N = len(ecc_bins)
    x = np.arange(N)

    def gaussian_fixed_mu(x, A, sigma, mu):
        return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

    params_fixed = np.zeros((N, 3))
    fitted_fixed = np.zeros_like(D)

    for i in range(N):
        y = D[i]

        # Force μ = i (diagonal constraint)
        mu = float(i)

        # Initial guesses
        A0 = y[i] if y[i] > 0 else np.max(y)
        sigma0 = 1.5

        try:
            popt, _ = curve_fit(
                lambda x, A, sigma: gaussian_fixed_mu(x, A, sigma, mu),
                x, y,
                p0=[A0, sigma0],
                bounds=([0, 1e-3], [np.inf, N]),   # sigma >= 0
            )
            A_fit, sigma_fit = popt

        except Exception:
            A_fit, sigma_fit = A0, sigma0

        params_fixed[i] = [A_fit, sigma_fit, mu]
        fitted_fixed[i] = gaussian_fixed_mu(x, A_fit, sigma_fit, mu)

    return params_fixed, fitted_fixed

def fit_matrix_rows_truncated(D, ecc_bins, fix_mu=False, normalize=False):
    """
    Fit truncated Gaussian to every row of D.

    Returns:
        params  (N×P)
        fitted  (N×N)
    """
    x = ecc_bin_centers(ecc_bins)
    n = D.shape[0]

    params = []
    fitted = np.zeros_like(D)

    for i in range(n):
        row = D[i]
        mu_init = x[i]

        popt, yfit = fit_row_truncated_gaussian(
            row, x, mu_init,
            fix_mu=fix_mu,
            normalize=normalize
        )
        params.append(popt)
        fitted[i] = yfit

    return np.array(params), fitted


# ============================================================
# === Plot Diagonal Observed vs Expected
# ============================================================

def plot_diag_expected_vs_observed(D, F, ecc_bins, out_png):
    """Scatter: expected vs observed on diagonal."""
    x = ecc_bin_centers(ecc_bins)
    obs = np.diag(D)
    exp = np.diag(F)

    plt.figure(figsize=(5,5))
    plt.scatter(exp, obs, c='k', alpha=0.7)
    mn = min(exp.min(), obs.min())
    mx = max(exp.max(), obs.max())
    plt.plot([mn, mx], [mn, mx], 'r--')
    plt.xlabel("Expected (Gaussian fit)")
    plt.ylabel("Observed (data)")
    plt.title("Diagonal: observed vs Gaussian expectation")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_diagonal_observed_vs_expected(D, F, ecc_labels, out_png):
    """
    Scatter plot comparing observed vs Gaussian-expected diagonal values.

    Parameters
    ----------
    D : np.ndarray
        Observed density matrix (NxN).
    F : np.ndarray
        Gaussian-fitted matrix (NxN).
    ecc_labels : list of str
        The eccentricity bin labels (same as used in plot_matrix).
    out_png : Path
        Output image path.

    What it shows:
    --------------
    - Each point = one eccentricity bin (diagonal element)
    - x-axis = expected value from the Gaussian model
    - y-axis = observed value in your data
    - Red dashed line = perfect match (y = x)

    Interpretation:
    ---------------
    Points near the line → model matches well  
    Points far from line → deviations from the Gaussian connectivity pattern
    """
    obs = np.diag(D)
    exp = np.diag(F)

    plt.figure(figsize=(6, 6))
    plt.scatter(exp, obs, c='k', alpha=0.8)

    # identity line
    mn = min(exp.min(), obs.min())
    mx = max(exp.max(), obs.max())
    plt.plot([mn, mx], [mn, mx], 'r--', label="y = x")

    # label each point with its ecc label
    for x, y, lbl in zip(exp, obs, ecc_labels):
        plt.text(x, y, lbl, fontsize=8, ha='right', va='bottom')

    plt.xlabel("Expected (Gaussian fit)")
    plt.ylabel("Observed (data)")
    plt.title("Diagonal Connectivity: Observed vs Gaussian Expected")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_eccentricity_profile_vs_model(D, F, ecc_labels, out_png):
    """
    Line plot showing observed and Gaussian-expected diagonal values 
    across eccentricity bins.

    Parameters
    ----------
    D : np.ndarray
        Observed density matrix.
    F : np.ndarray
        Gaussian-fitted matrix.
    ecc_labels : list of str
        Labels of eccentricity bins (same order as matrix rows/cols).
    out_png : Path
        Output file.

    What it shows:
    --------------
    - Blue line: observed diagonal connectivity for each eccentricity bin
    - Orange line: expected diagonal from Gaussian model
    - x-axis: eccentricity bins (same labels used in your heatmaps)

    Interpretation:
    ---------------
    This reveals how "self-connectivity" (bin→itself) varies as a 
    function of eccentricity, and whether the Gaussian fit respects 
    that dependence. Foveal vs peripheral differences appear clearly.
    """
    obs = np.diag(D)
    exp = np.diag(F)

    plt.figure(figsize=(7, 4))
    x = np.arange(len(ecc_labels))

    plt.plot(x, obs, '-o', label="Observed", linewidth=2)
    plt.plot(x, exp, '-o', label="Gaussian Expected", linewidth=2)

    plt.xticks(x, [e.replace("_", "–") for e in ecc_labels], rotation=45)
    plt.ylabel("Diagonal Connectivity (density)")
    plt.title("Eccentricity Profile: Observed vs Gaussian Model")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_gaussian_mu_and_diagonal_deviation(params, D, ecc_bins, out_dir):
    """
    params: array N×3 of [mu, sigma, A]
    D: observed matrix
    plots:
      1) expected mu_i vs eccentricity
      2) deviation Δ_i = observed_diag - expected_amp
    """

    mu = params[:,0]
    sigma = params[:,1]
    amp = params[:,2]

    # Observed diagonal
    diag_obs = np.diag(D)

    # Deviation from expected peak height
    deviation = diag_obs - amp

    x = np.arange(len(ecc_bins))
    labels = [b.replace("_", "–") for b in ecc_bins]

    # -----------------------------
    # Plot 1 — Gaussian μ vs eccentricity
    # -----------------------------
    plt.figure(figsize=(6,4))
    plt.plot(x, mu, 'o--', color='purple', label='Gaussian μ (peak location)')
    plt.plot(x, x, 'k:', label='Ideal diagonal (μ=i)')
    plt.xticks(x, labels, rotation=45)
    plt.xlabel("Eccentricity bin")
    plt.ylabel("Mean μ (bin index)")
    plt.title("Gaussian Peak Location per Eccentricity Bin")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "gaussian_mu_vs_eccentricity.png", dpi=300)
    plt.close()

    # -----------------------------
    # Plot 2 — Deviation of observed diagonal from expected Gaussian peak
    # -----------------------------
    plt.figure(figsize=(6,4))
    plt.axhline(0, color='gray', linestyle='--')
    plt.plot(x, deviation, 'o-', color='darkred', label='Deviation Δ = D[ii] - A_i')
    plt.xticks(x, labels, rotation=45)
    plt.xlabel("Eccentricity bin")
    plt.ylabel("Deviation from expected peak")
    plt.title("Observed vs Expected (Gaussian Peak Height)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "gaussian_diagonal_deviation.png", dpi=300)
    plt.close()

def plot_diag_observed_vs_expected_fixed_mu(D, F, ecc_bins, out_png):
    """
    Observed D[i,i] vs expected F[i,i] using forced-mu Gaussian fits.
    """

    obs = np.diag(D)
    exp = np.diag(F)

    plt.figure(figsize=(5,5))
    plt.scatter(exp, obs, c='purple', s=50, label="data")
    mn = min(exp.min(), obs.min())
    mx = max(exp.max(), obs.max())
    plt.plot([mn, mx], [mn, mx], 'k--')

    plt.xlabel("Expected (Gaussian, μ=i)")
    plt.ylabel("Observed (data)")
    plt.title("Observed vs expected diagonal (forced μ)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_gaussian_profile_fixed_mu(D, F, ecc_bins, out_png):
    """
    Plot for each row i:
       - data row D[i]
       - Gaussian fit F[i] centered at diagonal (μ=i)
    """

    x = np.arange(len(ecc_bins))
    plt.figure(figsize=(12, 6))

    for i in range(len(ecc_bins)):
        plt.plot(x, D[i], 'o-', alpha=0.5)
        plt.plot(x, F[i], '-', linewidth=2)

    plt.xticks(x, [e.replace("_", "–") for e in ecc_bins], rotation=45)
    plt.xlabel("Eccentricity")
    plt.ylabel("Density")
    plt.title("Row-wise Gaussian fits (μ fixed to diagonal)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_sigma_vs_eccentricity(params, ecc_labels, out_png):
    """
    Plot sigma (spread of the Gaussian) for each row.
    params: Nx3 array from fit_matrix_rows → columns = [A, mu, sigma]
    """
    import matplotlib.pyplot as plt
    import numpy as np

    sigmas = params[:, 2]  # column 2 = sigma

    x = np.arange(len(sigmas))

    plt.figure(figsize=(8,4))
    plt.plot(x, sigmas, 'o-', color='purple')

    plt.xticks(x, ecc_labels, rotation=45)
    plt.ylabel("Sigma (spread)")
    plt.xlabel("Eccentricity bin")
    plt.title("Gaussian σ across eccentricity bins")

    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_deviation_from_expected(diagonal_obs, diagonal_exp, sigma, ecc_labels, out_png):
    """
    Plot deviation Δ(i) = observed[i] - expected[i]
    plus Gaussian sigma as vertical error bars.

    sigma[i] = fitted Gaussian width for row i.
    """

    import matplotlib.pyplot as plt
    import numpy as np

    delta = diagonal_obs - diagonal_exp
    x = np.arange(len(delta))

    plt.figure(figsize=(8, 4))

    # Zero line = perfect match
    plt.axhline(0, color='gray', linestyle='--', linewidth=1)

    # Error bars for sigma
    plt.errorbar(
        x, delta, yerr=sigma,
        fmt='o', color='darkred', ecolor='black', elinewidth=1.2,
        capsize=4, label="Deviation Δ with σ error"
    )

    plt.plot(x, delta, color='darkred', linewidth=1)

    plt.xticks(x, ecc_labels, rotation=45)
    plt.ylabel("Deviation from Gaussian peak (Observed − Expected)")
    plt.title("Diagonal deviation from expected Gaussian peak")

    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()



def plot_deviation_from_expected_fixed_mu(diagonal_obs, diagonal_exp, sigma, ecc_labels, out_png):
    """
    Deviation Δ(i) = observed[i] - expected_fixed_mu[i]
    using the FIXED-MU Gaussian fits (mu = diagonal index),
    with sigma displayed as vertical error bars.
    """

    import matplotlib.pyplot as plt
    import numpy as np

    delta = diagonal_obs - diagonal_exp
    x = np.arange(len(delta))

    plt.figure(figsize=(8, 4))

    # Zero line
    plt.axhline(0, color='gray', linestyle='--', linewidth=1)

    # Sigma error bars
    plt.errorbar(
        x, delta, yerr=sigma,
        fmt='o', color='darkblue', ecolor='black', elinewidth=1.2,
        capsize=4, label="Deviation Δ with σ error"
    )

    plt.plot(x, delta, color='darkblue', linewidth=1)

    plt.xticks(x, ecc_labels, rotation=45)
    plt.ylabel("Deviation (Observed − Expected fixed-μ)")
    plt.title("Diagonal deviation from forced-on-diagonal Gaussian peak (with σ)")

    plt.grid(alpha=0.25)
    plt.tight_layout()
    save_figure(out_png, dpi=300)
    plt.close()


def plot_diag_observed_vs_expected_truncated(D, F, ecc_labels, out_png):
    """
    Observed vs expected diagonal for truncated Gaussian models.
    - D: observed matrix
    - F: fitted truncated Gaussian matrix
    """

    import numpy as np
    import matplotlib.pyplot as plt

    obs = np.diag(D)
    exp = np.diag(F)

    plt.figure(figsize=(5,5))
    plt.scatter(exp, obs, c="k", s=60, label="Bins")
    mn = min(exp.min(), obs.min())
    mx = max(exp.max(), obs.max())
    plt.plot([mn, mx], [mn, mx], "r--", label="y = x")

    plt.xlabel("Expected (Truncated Gaussian)")
    plt.ylabel("Observed")
    plt.title("Diagonal: observed vs truncated-Gaussian prediction")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    save_figure(out_png, dpi=300)
    plt.close()

def plot_gaussian_profile_truncated(D, F, ecc_labels, out_png):
    """
    Plot each row's observed profile and its truncated-Gaussian fit.
    """

    import numpy as np
    import matplotlib.pyplot as plt

    x = np.arange(len(ecc_labels))

    plt.figure(figsize=(10, 6))

    for i in range(D.shape[0]):
        plt.plot(x, D[i], "k-", alpha=0.3)
        plt.plot(x, F[i], "-", alpha=0.8)

    plt.xticks(x, ecc_labels, rotation=45)
    plt.ylabel("Value")
    plt.title("Observed vs truncated-Gaussian fit (row profiles)")
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_deviation_from_expected_truncated(
    diagonal_obs,
    diagonal_exp,
    sigma,
    ecc_labels,
    out_png,
    xlim=None   # <-- NEW optional axis limits (e.g., (-10, 10))
):
    """
    Horizontal deviation plot for truncated Gaussian with FREE μ.
    Δ = observed - expected.
    """


    delta = diagonal_obs - diagonal_exp
    y = np.arange(len(delta))

    # Pretty labels
    nice_labels = [lbl.replace("_", "–") + "°" for lbl in ecc_labels]

    fig, ax = plt.subplots(figsize=(4, 8))

    # Clean style: remove top/right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Vertical zero line
    ax.axvline(0, color="lightgray", linewidth=2, zorder=1)

    # Error bars – horizontal, no caps
    ax.errorbar(
        delta,                    # x-values
        y,                        # y-values
        xerr=sigma,               # horizontal σ
        fmt="o",
        markersize=6,
        color="darkred",
        ecolor="gray",
        elinewidth=2,
        capsize=0,                # <-- NO T-shaped ends
        zorder=3,
    )

    ax.set_yticks(y)
    ax.set_yticklabels(nice_labels)
    ax.set_xlabel("Deviation (Observed – Expected)")
    ax.set_title("Diagonal deviation — Truncated Gaussian (μ free)")

    # Optional X limits
    if xlim is not None:
        ax.set_xlim(xlim)

    plt.tight_layout()
    save_figure(out_png, dpi=300)
    plt.close()




def plot_diag_observed_vs_expected_truncated_fixed_mu(D, F, ecc_labels, out_png):
    """
    Observed vs expected on the diagonal for truncated Gaussian
    with mu forced to i (the diagonal index).
    """

    obs = np.diag(D)
    exp = np.diag(F)

    plt.figure(figsize=(5,5))
    plt.scatter(exp, obs, c="navy", s=60)
    mn = min(exp.min(), obs.min())
    mx = max(exp.max(), obs.max())
    plt.plot([mn, mx], [mn, mx], "r--")

    plt.xlabel("Expected (Truncated Gaussian, fixed μ=i)")
    plt.ylabel("Observed")
    plt.title("Diagonal: observed vs expected (truncated, fixed μ)")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    save_figure(out_png, dpi=300)
    plt.close()

def plot_gaussian_profile_truncated_fixed_mu(D, F, ecc_labels, out_png):
    """
    Plot each row's observed profile and the truncated-Gaussian fit
    with μ=i enforced.
    """
    x = np.arange(len(ecc_labels))

    plt.figure(figsize=(10, 6))

    for i in range(D.shape[0]):
        plt.plot(x, D[i], "k-", alpha=0.3)
        plt.plot(x, F[i], "-", color="navy", alpha=0.8)

    plt.xticks(x, ecc_labels, rotation=45)
    plt.ylabel("Value")
    plt.title("Observed vs truncated Gaussian (fixed μ=i) row profiles")
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_deviation_from_expected_truncated_fixed_mu(
    diagonal_obs,
    diagonal_exp,
    sigma,
    ecc_labels,
    out_png,
    xlim=None   # <-- NEW optional axis limits
):
    """
    Horizontal deviation plot for truncated Gaussian with FIXED μ.
    Δ = observed - expected_fixed.
    """

    import numpy as np
    import matplotlib.pyplot as plt

    delta = diagonal_obs - diagonal_exp
    y = np.arange(len(delta))

    # Pretty labels
    nice_labels = [lbl.replace("_", "–") + "°" for lbl in ecc_labels]

    fig, ax = plt.subplots(figsize=(4, 8))

    # Clean style
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Zero deviation line
    ax.axvline(0, color="lightgray", linewidth=2, zorder=1)

    # Horizontal error bars
    ax.errorbar(
        delta,
        y,
        xerr=sigma,
        fmt="o",
        markersize=6,
        color="navy",
        ecolor="gray",
        elinewidth=2,
        capsize=0,
        zorder=3,
    )

    ax.set_yticks(y)
    ax.set_yticklabels(nice_labels)
    ax.set_xlabel("Deviation (Observed – Expected, μ fixed)")
    ax.set_title("Diagonal deviation — Truncated Gaussian (μ fixed)")

    # Optional limits
    if xlim is not None:
        ax.set_xlim(xlim)

    plt.tight_layout()
    save_figure(out_png, dpi=300)
    plt.close()

def plot_mean_deviation_with_sem(
    diagonal_obs,
    diagonal_exp,
    out_png,
    title="Mean deviation from model (Observed − Expected)"
):
    """
    Plot the MEAN deviation and SEM across all eccentricity bins.
    Produces a single point with vertical error bar.

    Works with:
      - truncated Gaussian
      - non-truncated Gaussian
      - μ fixed or free
    """

    import numpy as np
    import matplotlib.pyplot as plt

    delta = diagonal_obs - diagonal_exp
    mean_delta = np.mean(delta)
    sem_delta  = np.std(delta, ddof=1) / np.sqrt(len(delta))

    fig, ax = plt.subplots(figsize=(4, 6))

    # zero-line (perfect agreement)
    ax.axhline(0, color="lightgray", linewidth=2)

    ax.errorbar(
        x=[0],
        y=[mean_delta],
        yerr=[sem_delta],
        fmt="o",
        color="black",
        ecolor="gray",
        elinewidth=2,
        capsize=0,  # no caps = no "T"
        markersize=8,
    )

    ax.set_xticks([])
    ax.set_ylabel("Deviation (Observed − Expected)")
    ax.set_title(title)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    save_figure(out_png, dpi=300)
    plt.close()



def plot_sigma_vs_eccentricity_truncated_fixed_mu(params, ecc_labels, out_png):
    """
    Plot σ (width) per eccentricity for truncated Gaussian (fixed μ).
    Clean scientific style: no grid, no top/right frame.
    """
    sigma = params[:, 1]  # params = [A, σ, c] for fixed μ

    x = np.arange(len(sigma))

    # Nice eccentricity labels
    nice_labels = [lbl.replace("_", "–") + "°" for lbl in ecc_labels]

    plt.figure(figsize=(8,4))

    # Remove top/right borders
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Plot sigma
    plt.plot(x, sigma, 'o-', color='darkgreen', linewidth=1.8)

    plt.xticks(x, nice_labels, rotation=45)
    plt.ylabel("σ (spread)")
    plt.title("σ across eccentricity (Truncated Gaussian, fixed μ)")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_sigma_vs_eccentricity_truncated(params, ecc_labels, out_png):
    """
    Plot σ (width) per eccentricity for truncated Gaussian (free μ).
    Clean scientific style: no grid, no top/right frame.
    """

    sigma = params[:, 2]  # params = [A, μ, σ, c]

    x = np.arange(len(sigma))

    # Format labels: "0-2°"
    nice_labels = [lbl.replace("_", "–") + "°" for lbl in ecc_labels]

    plt.figure(figsize=(8,4))

    # Remove top/right borders
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Plot sigma
    plt.plot(x, sigma, 'o-', color='purple', linewidth=1.8)

    plt.xticks(x, nice_labels, rotation=45)
    plt.ylabel("σ (spread)")
    plt.title("σ across eccentricity (Truncated Gaussian, free μ)")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


# ============================================================
# === CLI Interface
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Fit Gaussian models to each row of a connectivity matrix.")
    parser.add_argument("--matrix", type=str, required=True, help="CSV matrix file (NxN).")
    parser.add_argument("--ecc-bins", type=str, required=True,
                        help="Comma-separated eccentricity bins, e.g. 0_1,1_2,2_3,...")
    parser.add_argument("--fix-mu", action="store_true",
                        help="Fix Gaussian center μ at the diagonal.")
    parser.add_argument("--out-dir", type=str, default="gaussian_fit_output",
                        help="Output directory.")

    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ecc_bins = [b.strip() for b in args.ecc_bins.split(",")]

    # Load matrix
    D = np.loadtxt(args.matrix, delimiter=",")

    # Fit
    params, fitted = fit_matrix_rows(D, ecc_bins, fix_mu=args.fix_mu)

    # Save
    np.savetxt(out_dir / "gaussian_params.csv", params, delimiter=",",
               header="k,mu,sigma,c" if not args.fix_mu else "k,sigma,c",
               comments='')

    np.savetxt(out_dir / "gaussian_fitted_matrix.csv", fitted, delimiter=",")

    # Plot comparison on diagonal
    plot_diag_expected_vs_observed(D, fitted, ecc_bins, out_dir / "diag_observed_vs_expected.png")

    print("\nGaussian fitting completed.")
    print(f"Saved parameters, reconstructed matrix, and plot to: {out_dir}\n")


if __name__ == "__main__":
    main()
