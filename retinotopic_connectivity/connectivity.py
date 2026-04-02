import subprocess
import numpy as np
import nibabel as nib
from pathlib import Path
from typing import List, Optional
from filelock import FileLock
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.colors as mcolors
# vendored modules (placed at repo root)
from utils import save_figure
from fit_gaussian_connectivity import (
    fit_matrix_rows,
    plot_diagonal_observed_vs_expected,
    plot_eccentricity_profile_vs_model,
    plot_gaussian_mu_and_diagonal_deviation,
    fit_matrix_rows_fixed_mu,
    plot_diag_observed_vs_expected_fixed_mu,
    plot_gaussian_profile_fixed_mu,
    plot_deviation_from_expected,
    plot_deviation_from_expected_fixed_mu,
    plot_sigma_vs_eccentricity,
    fit_matrix_rows_truncated,
    plot_diag_observed_vs_expected_truncated,
    plot_gaussian_profile_truncated,
    plot_deviation_from_expected_truncated,
    plot_diag_observed_vs_expected_truncated_fixed_mu,
    plot_gaussian_profile_truncated_fixed_mu,
    plot_deviation_from_expected_truncated_fixed_mu,
    plot_sigma_vs_eccentricity_truncated_fixed_mu,
    plot_sigma_vs_eccentricity_truncated,
    plot_mean_deviation_with_sem,
)
from diagonal_to_braplot import compute_shell_vals, plot_radial_shells
from matplotlib import cm, colors as mpl_colors
from matrix_to_dva import plot_dva_bar

# ---------------------------------------------------------------------
# Visual area labels (Benson-style)
# ---------------------------------------------------------------------
AREA_LABELS = ["V1", "V2", "V3", "hV4", "VO1", "VO2", "LO1", "LO2", "TO1", "TO2", "V3b", "V3a"]


def format_colorbar(cbar, data, vmin, vmax, log_scale):
    if log_scale:
        vmin_eff = vmin if vmin is not None else float(np.nanmin(data))
        vmax_eff = vmax if vmax is not None else float(np.nanmax(data))
        mid = np.sqrt(vmin_eff * vmax_eff) if vmin_eff > 0 else (vmax_eff / 10.0)
    else:
        vmin_eff = vmin if vmin is not None else float(np.nanmin(data))
        vmax_eff = vmax if vmax is not None else float(np.nanmax(data))
        mid = (vmin_eff + vmax_eff) / 2.0

    ticks = [vmin_eff, mid, vmax_eff]
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter("%.2g"))
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label("Streamline density", fontsize=12)


def load_mask(path: Path):
    return np.squeeze(nib.load(str(path)).get_fdata())


_area_cache = {}


def mask_area(mask_a: Path, mask_b: Path) -> float:
    """
    Area/volume of intersection between two binary masks, in mm^3.
    Uses voxel volume from mask_a header.
    """
    key = f"{mask_a.name}_{mask_b.name}"
    if key in _area_cache:
        return _area_cache[key]

    if (not mask_a.exists()) or (not mask_b.exists()):
        _area_cache[key] = 0.0
        return 0.0

    try:
        a = load_mask(mask_a) > 0
        b = load_mask(mask_b) > 0

        min_shape = tuple(np.minimum(a.shape, b.shape))
        slicer = tuple(slice(0, m) for m in min_shape)
        a = a[slicer]
        b = b[slicer]
        inter = a & b

        voxel_vol = np.prod(nib.load(str(mask_a)).header.get_zooms()[:3])
        area = float(np.sum(inter) * voxel_vol)
    except Exception:
        area = 0.0

    _area_cache[key] = area
    return area


def intersect_masks(mask1: Path, mask2: Path, out_path: Path) -> Path:
    """
    Intersection = mask1 & mask2.
    Thread/process safe via file lock.
    """
    lock = FileLock(str(out_path) + ".lock")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with lock:
        if out_path.exists():
            return out_path

        m1_img = nib.load(str(mask1))
        m2_img = nib.load(str(mask2))

        m1 = np.squeeze(m1_img.get_fdata()) > 0
        m2 = np.squeeze(m2_img.get_fdata()) > 0

        if len(m1.shape) != 3 or len(m2.shape) != 3:
            raise ValueError(f"Mask shapes must be 3D after squeeze: {m1.shape} vs {m2.shape}")

        min_shape = tuple(np.minimum(m1.shape, m2.shape))
        slicer = tuple(slice(0, m) for m in min_shape)
        inter = (m1[slicer] & m2[slicer]).astype(np.uint8)

        out_img = nib.Nifti1Image(inter, m1_img.affine, m1_img.header)
        out_img.set_data_dtype(np.uint8)
        nib.save(out_img, str(out_path))

    return out_path


def subject_threshold_map(base_map: Path, low: float, high: float, var_type=None):
    img = nib.load(str(base_map))
    data = img.get_fdata()
    if var_type == "angle":
        data = np.abs(data)
    mask = (data >= low) & (data <= high)
    return mask.astype(np.uint8), img


def make_subject_patch_mask(
    ecc_map: Path,
    ang_map: Optional[Path],
    ecc_range: str,
    ang_range: Optional[str],
    out_dir: Path,
) -> Path:
    """
    Create an ROI mask for a subject map:
      - ecc-only if ang_range is None
      - ecc x polar if ang_range provided

    Output is saved under out_dir.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    if ang_range is None:
        out = out_dir / f"ecc_{ecc_range}.nii.gz"
        if out.exists():
            return out
        mask, img = subject_threshold_map(ecc_map, *map(float, ecc_range.split("_")))
        nib.save(nib.Nifti1Image(mask, img.affine, img.header), str(out))
        return out

    out = out_dir / f"ecc_{ecc_range}_polar_{ang_range}.nii.gz"
    if out.exists():
        return out

    ecc_low, ecc_high = map(float, ecc_range.split("_"))
    ang_low, ang_high = map(float, ang_range.split("_"))

    ecc_mask, img = subject_threshold_map(ecc_map, ecc_low, ecc_high)
    ang_mask, _ = subject_threshold_map(ang_map, ang_low, ang_high, var_type="angle")

    roi = ((ecc_mask > 0) & (ang_mask > 0)).astype(np.uint8)
    nib.save(nib.Nifti1Image(roi, img.affine, img.header), str(out))
    return out


def extract_visual_area_mask(varea_img: Path, area_name: str, out_path: Path) -> Path:
    """
    Extract a binary mask for a visual area from varea.nii.gz assuming Benson labels.
    """
    if out_path.exists():
        return out_path

    if area_name not in AREA_LABELS:
        raise ValueError(f"Unknown area '{area_name}'. Valid: {AREA_LABELS}")

    val = AREA_LABELS.index(area_name) + 1  # 1-based

    img = nib.load(str(varea_img))
    data = img.get_fdata()
    mask = (data == val).astype(np.uint8)
    nib.save(nib.Nifti1Image(mask, img.affine, img.header), str(out_path))
    return out_path


def run_tckedit(track: Path, roi1: Path, roi2: Path, out_tck: Path, ends_only: bool = True, roi_order: bool = False) -> int:
    out_tck.parent.mkdir(parents=True, exist_ok=True)

    if not roi_order:
        cmd = ["tckedit", str(track), str(out_tck), "-include", str(roi1), "-include", str(roi2)]
        if ends_only:
            cmd.append("-ends_only")
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    else:
        tmp12 = out_tck.with_suffix(".tmp12.tck")
        tmp21 = out_tck.with_suffix(".tmp21.tck")

        cmd12 = ["tckedit", str(track), str(tmp12), "-include_ordered", str(roi1), "-include_ordered", str(roi2)]
        cmd21 = ["tckedit", str(track), str(tmp21), "-include_ordered", str(roi2), "-include_ordered", str(roi1)]
        if ends_only:
            cmd12.append("-ends_only")
            cmd21.append("-ends_only")

        subprocess.run(cmd12, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(cmd21, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(["tckedit", str(tmp12), str(tmp21), str(out_tck)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        tmp12.unlink(missing_ok=True)
        tmp21.unlink(missing_ok=True)

    try:
        if not out_tck.exists():
            return 0
        c = int(subprocess.check_output(["tckinfo", str(out_tck), "-count"]).decode().split()[-1])
    except Exception:
        c = 0
    return c


def plot_matrix(matrix, title, path: Path, labels, cmap="viridis", vmin=None, vmax=None, log_scale=False):
    plt.figure(figsize=(5, 4))
    mat = matrix.astype(float)

    if log_scale:
        mat = np.where(mat <= 0, 1e-6, mat)
        im = plt.imshow(mat, cmap=cmap, norm=LogNorm(vmin=vmin or float(np.nanmin(mat)), vmax=vmax or float(np.nanmax(mat))))
    else:
        im = plt.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(im, label="Streamline density")
    format_colorbar(cbar, mat, vmin, vmax, log_scale)

    plt.xticks(range(len(labels)), [e.replace("_", "–") for e in labels], rotation=45)
    plt.yticks(range(len(labels)), [e.replace("_", "–") for e in labels])
    plt.xlabel("Eccentricity bins (DVA)")
    plt.ylabel("Eccentricity bins (DVA)")
    plt.title(title + (" (Log Scale)" if log_scale else ""), fontsize=10)
    plt.tight_layout()
    save_figure(path, dpi=300)
    plt.close()


def plot_dva_bar_save(shell_vals, out_png: Path, bar_x_dim: float = 100, base_width: float = 6, base_height: float = 4):
    labels = ["same DVA", "±1 DVA", "±2 DVA", "±3 DVA", "±4 DVA", "±5 DVA", "≥6 DVA"]

    means = [shell_vals[k].mean() for k in range(6)]
    sems  = [shell_vals[k].std() / np.sqrt(len(shell_vals[k])) for k in range(6)]
    far_vals = np.concatenate([v for k, v in shell_vals.items() if k >= 6])
    means.append(far_vals.mean())
    sems.append(far_vals.std() / np.sqrt(len(far_vals)))

    n_bars = len(means)
    shrink = n_bars / bar_x_dim
    fig_w = base_width * shrink if shrink < 1 else base_width
    fig_h = base_height

    bar_w = 0.70 * shrink if shrink < 1 else 0.70
    x = np.arange(n_bars) * shrink

    base_cmap = cm.hot
    trunc_cmap = mpl_colors.LinearSegmentedColormap.from_list("hot_truncated", base_cmap(np.linspace(0.0, 0.85, 256)))
    norm = mpl_colors.Normalize(vmin=min(means), vmax=max(means))
    bar_colors = list(trunc_cmap(norm(means)))

    def lighten(color, amount=0.35):
        r, g, b, a = color
        return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount, a)

    bar_colors[0] = lighten(bar_colors[0])
    bar_colors[1] = lighten(bar_colors[1])

    plt.figure(figsize=(fig_w, fig_h))
    plt.bar(x, means, yerr=sems, width=bar_w, color=bar_colors, capsize=0, ecolor="0.5", error_kw=dict(linewidth=1))

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ymin = 0
    ymax = max(np.array(means) + np.array(sems))
    ax.set_ylim(ymin, ymax * 1.05)
    ax.set_yticks([ymin, ymax / 2, ymax])
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30)
    plt.ylabel("Mean streamline density")

    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close()


def run_all_gaussian_fitting(D, ecc_bins, out_dir: Path, fit_truncated_gaussian_normalized: bool):
    out_dir.mkdir(parents=True, exist_ok=True)

    # Gaussian free mu
    params, fitted = fit_matrix_rows(D, ecc_bins, fix_mu=False)
    gdir = out_dir / "gaussian_fit"
    gdir.mkdir(exist_ok=True)
    np.savetxt(gdir / "gaussian_params.csv", params, delimiter=",")
    np.savetxt(gdir / "gaussian_fitted_matrix.csv", fitted, delimiter=",")
    plot_diagonal_observed_vs_expected(D, fitted, ecc_bins, gdir / "diag_observed_vs_expected.png")
    plot_eccentricity_profile_vs_model(D, fitted, ecc_bins, gdir / "eccentricity_profile_obs_vs_model.png")
    plot_gaussian_mu_and_diagonal_deviation(params, D, ecc_bins, gdir)

    diag_obs = np.diag(D)
    diag_free = np.diag(fitted)
    sigma_free = params[:, 2]
    plot_deviation_from_expected(diag_obs, diag_free, sigma_free, ecc_bins, gdir / "deviation_centerline_free_mu.png")
    plot_sigma_vs_eccentricity(params, ecc_bins, gdir / "sigma_free_mu.png")

    # Gaussian fixed mu
    params_fixed, fitted_fixed = fit_matrix_rows_fixed_mu(D, ecc_bins)
    gdir = out_dir / "gaussian_fit_fixed_mu"
    gdir.mkdir(exist_ok=True)
    np.savetxt(gdir / "params_fixed_mu.csv", params_fixed, delimiter=",")
    np.savetxt(gdir / "fitted_matrix_fixed_mu.csv", fitted_fixed, delimiter=",")
    plot_diag_observed_vs_expected_fixed_mu(D, fitted_fixed, ecc_bins, gdir / "diag_observed_vs_expected_fixed_mu.png")
    plot_gaussian_profile_fixed_mu(D, fitted_fixed, ecc_bins, gdir / "profile_fixed_mu.png")
    A_fixed = params_fixed[:, 0]
    sigma_fixed = params_fixed[:, 1]
    plot_deviation_from_expected_fixed_mu(np.diag(D), A_fixed, sigma_fixed, ecc_bins, gdir / "deviation_fixed_mu_with_sigma.png")
    plot_sigma_vs_eccentricity(params_fixed, ecc_bins, gdir / "sigma_fixed_mu.png")

    # Truncated Gaussian free mu
    gdir = out_dir / "gaussian_truncated_fit"
    gdir.mkdir(exist_ok=True)
    params_trunc_free, fitted_trunc_free = fit_matrix_rows_truncated(D, ecc_bins, fix_mu=False, normalize=fit_truncated_gaussian_normalized)
    np.savetxt(gdir / "params_truncated_free_mu.csv", params_trunc_free, delimiter=",")
    np.savetxt(gdir / "fitted_matrix_truncated_free_mu.csv", fitted_trunc_free, delimiter=",")
    diag_obs = np.diag(D)
    diag_exp = np.diag(fitted_trunc_free)
    sigma_free = params_trunc_free[:, 2]
    plot_diag_observed_vs_expected_truncated(D, fitted_trunc_free, ecc_bins, gdir / "diag_observed_vs_expected_truncated_free_mu.png")
    plot_gaussian_profile_truncated(D, fitted_trunc_free, ecc_bins, gdir / "profile_truncated_free_mu.png")
    plot_deviation_from_expected_truncated(diag_obs, diag_exp, sigma_free, ecc_bins, gdir / "deviation_truncated_free_mu_with_sigma.png", xlim=(-10, 10))
    plot_sigma_vs_eccentricity_truncated(params_trunc_free, ecc_bins, gdir / "sigma_truncated_free_mu.png")
    plot_mean_deviation_with_sem(diag_obs, diag_exp, gdir / "mean_deviation_truncated_free_mu.png", title="Mean deviation — Truncated Gaussian (μ free)")

    # Truncated Gaussian fixed mu
    gdir = out_dir / "gaussian_truncated_fit_fixed_mu"
    gdir.mkdir(exist_ok=True)
    params_trunc_fixed, fitted_trunc_fixed = fit_matrix_rows_truncated(D, ecc_bins, fix_mu=True, normalize=fit_truncated_gaussian_normalized)
    np.savetxt(gdir / "params_truncated_fixed_mu.csv", params_trunc_fixed, delimiter=",")
    np.savetxt(gdir / "fitted_matrix_truncated_fixed_mu.csv", fitted_trunc_fixed, delimiter=",")
    diag_fixed = np.diag(fitted_trunc_fixed)
    sigma_fixed = params_trunc_fixed[:, 1]
    plot_diag_observed_vs_expected_truncated_fixed_mu(D, fitted_trunc_fixed, ecc_bins, gdir / "diag_observed_vs_expected_truncated_fixed_mu.png")
    plot_gaussian_profile_truncated_fixed_mu(D, fitted_trunc_fixed, ecc_bins, gdir / "profile_truncated_fixed_mu.png")
    plot_deviation_from_expected_truncated_fixed_mu(np.diag(D), diag_fixed, sigma_fixed, ecc_bins, gdir / "deviation_truncated_fixed_mu_with_sigma.png", xlim=(-10, 10))
    plot_sigma_vs_eccentricity_truncated_fixed_mu(params_trunc_fixed, ecc_bins, gdir / "sigma_truncated_fixed_mu.png")
    plot_mean_deviation_with_sem(diag_obs, diag_fixed, gdir / "mean_deviation_truncated_fixed_mu.png", title="Mean deviation — Truncated Gaussian (μ fixed)")


def run_tckedit_endpoints_in_mask(track: Path, mask: Path, out_tck: Path, ends_only: bool = True) -> None:
    """
    Filter *track* so that both streamline endpoints lie within *mask*.
    Equivalent to: tckedit track out_tck -include mask -include mask [-ends_only]
    """
    out_tck.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["tckedit", str(track), str(out_tck), "-include", str(mask), "-include", str(mask), "-force"]
    if ends_only:
        cmd.append("-ends_only")
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)


def run_tck2connectome(tck: Path, labels: Path, out_csv: Path, n_nodes: int = 12) -> None:
    """
    Run MRtrix tck2connectome to produce a *n_nodes* × *n_nodes* area connectivity CSV.
    Flags match original VISCONTI script:
      -symmetric -zero_diagonal -assignment_end_voxels -scale_invnodevol -force
    """
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "tck2connectome",
        str(tck),
        str(labels),
        str(out_csv),
        "-symmetric",
        "-zero_diagonal",
        "-assignment_end_voxels",
        "-scale_invnodevol",
        "-force",
    ]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)


def make_masked_label_image(varea_img: Path, mask: Path, out_path: Path) -> Path:
    """
    Multiply *varea_img* by *mask* to produce a label image restricted to the mask.
    Preserves affine and header from *varea_img*.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lock = FileLock(str(out_path) + ".lock")
    with lock:
        if out_path.exists():
            return out_path
        varea_nib = nib.load(str(varea_img))
        mask_data = np.squeeze(nib.load(str(mask)).get_fdata()) > 0
        varea_data = np.squeeze(varea_nib.get_fdata())
        masked = (varea_data * mask_data.astype(varea_data.dtype))
        out_img = nib.Nifti1Image(masked, varea_nib.affine, varea_nib.header)
        nib.save(out_img, str(out_path))
    return out_path


def plot_area_matrix(matrix: np.ndarray, title: str, path: Path, labels: List[str], cmap: str = "viridis", vmin: Optional[float] = None, vmax: Optional[float] = None, log_scale: bool = False):
    """Plot a 12×12 visual-area connectivity matrix."""
    plt.figure(figsize=(7, 6))
    mat = matrix.astype(float)

    if log_scale:
        mat = np.where(mat <= 0, 1e-6, mat)
        im = plt.imshow(mat, cmap=cmap, norm=LogNorm(vmin=vmin or float(np.nanmin(mat)), vmax=vmax or float(np.nanmax(mat))))
    else:
        im = plt.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(im, label="Connectivity (invnodevol)")
    format_colorbar(cbar, mat, vmin, vmax, log_scale)

    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.yticks(range(len(labels)), labels)
    plt.xlabel("Visual area")
    plt.ylabel("Visual area")
    plt.title(title + (" (Log Scale)" if log_scale else ""), fontsize=10)
    plt.tight_layout()
    save_figure(path, dpi=300)
    plt.close()
    
def make_smooth_colormap(color, name="custom", n_colors=256, light_mix=0.9, dark_mix=0.1):
    """
    Create a smooth colormap around a base color.

    Parameters
    ----------
    color : str | hex | RGB tuple
        Base color.
    name : str
        Name of the colormap.
    n_colors : int
        Number of gradient steps.
    light_mix : float
        How much white to mix for the light end.
    dark_mix : float
        How much black to mix for the dark end.
    """

    base = np.array(mcolors.to_rgb(color))
    white = np.array([1, 1, 1])
    black = np.array([0, 0, 0])

    light = base * (1 - light_mix) + white * light_mix
    dark = base * (1 - dark_mix) + black * dark_mix

    colors = [light, base, dark]

    cmap = mcolors.LinearSegmentedColormap.from_list(name, colors, N=n_colors)

    return cmap

def run_areas_per_bin_connectome(
    tract_tck: Path,
    ecc_map: Path,
    polar_map: Optional[Path],
    varea_map: Path,
    ecc_bins: List[str],
    polar_bins: List[str],
    outdir: Path,
    ends_only: bool,
    color_map: str,
    log_scale: bool,
    vmin: Optional[float],
    vmax: Optional[float],
):
    """
    Areas-per-eccentricity mode using MRtrix tck2connectome.

    For each eccentricity bin (and optionally each polar bin):
      1. create a combined retinotopic mask
      2. restrict the visual-area label image to that mask
      3. filter streamlines so both endpoints lie within the mask
      4. run tck2connectome -assignment_end_voxels

    Outputs are stored under: outdir/areas_per_bin_connectome/
    """
    n_areas = len(AREA_LABELS)
    ape_dir = outdir / "areas_per_bin_connectome"
    ape_dir.mkdir(parents=True, exist_ok=True)

    roi_dir = ape_dir / "ROIs"
    roi_dir.mkdir(exist_ok=True)

    tck_dir = ape_dir / "tcks"
    tck_dir.mkdir(exist_ok=True)

    color_map_list = [
        make_smooth_colormap(c.strip(), name=f"ecc_smooth_{i}")
        for i, c in enumerate(color_map.split(","))
    ]

    for idx, ecc in enumerate(ecc_bins):
        cmap = color_map_list[idx]

        for polar in polar_bins:
            use_polar = polar.lower() != "all"

            # --- combined mask (ecc only or ecc x polar) ---
            if use_polar:
                combined_mask = make_subject_patch_mask(
                    ecc_map, polar_map, ecc, polar, roi_dir / "ecc_polar"
                )
                tag = f"ecc{ecc}_polar{polar}"
            else:
                combined_mask = make_subject_patch_mask(
                    ecc_map, None, ecc, None, roi_dir / "ecc"
                )
                tag = f"ecc{ecc}"

            # --- masked label image (varea restricted to combined mask) ---
            label_img = make_masked_label_image(
                varea_map, combined_mask, roi_dir / f"labels_{tag}.nii.gz"
            )

            # --- filter tractogram: both endpoints inside combined mask ---
            filtered_tck = tck_dir / f"filtered_{tag}.tck"
            if not filtered_tck.exists():
                run_tckedit_endpoints_in_mask(
                    tract_tck, combined_mask, filtered_tck, ends_only=ends_only
                )

            # --- tck2connectome ---
            out_csv = ape_dir / f"area_matrix_{tag}.csv"
            if filtered_tck.exists() and not out_csv.exists():
                run_tck2connectome(filtered_tck, label_img, out_csv, n_nodes=n_areas)

            # --- load matrix and plot ---
            try:
                mat = np.loadtxt(str(out_csv), delimiter=",")
                if mat.shape != (n_areas, n_areas):
                    mat = np.zeros((n_areas, n_areas))
            except Exception:
                mat = np.zeros((n_areas, n_areas))

            ecc_label = ecc.replace("_", "–")
            if use_polar:
                polar_label = polar.replace("_", "–")
                title = f"Area connectivity  ecc {ecc_label}°  polar {polar_label}°"
            else:
                title = f"Area connectivity  ecc {ecc_label}°"

            plot_area_matrix(
                mat,
                title,
                ape_dir / f"area_matrix_{tag}",
                AREA_LABELS,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                log_scale=log_scale,
            )


def run_areas_per_bin_pairwise(
    tract_tck: Path,
    ecc_map: Path,
    polar_map: Optional[Path],
    varea_map: Path,
    ecc_bins: List[str],
    polar_bins: List[str],
    outdir: Path,
    ends_only: bool,
    roi_order: bool,
    color_map: str,
    log_scale: bool,
    vmin: Optional[float],
    vmax: Optional[float],
    zero_diagonal: bool = True,
    symmetric = True,
):
    """
    Areas-per-eccentricity mode using explicit pairwise ROI counting via run_tckedit().

    For each eccentricity bin (and optionally each polar bin):
      1. create a combined retinotopic mask
      2. intersect that mask with each visual area
      3. explicitly count streamlines between every pair of visual areas
         using run_tckedit(..., ends_only=..., roi_order=...)
      4. build the 12x12 density matrix manually

    Outputs are stored under: outdir/areas_per_bin_pairwise/
    """
    n_areas = len(AREA_LABELS)
    ape_dir = outdir / "areas_per_bin_pairwise"
    ape_dir.mkdir(parents=True, exist_ok=True)

    roi_dir = ape_dir / "ROIs"
    roi_dir.mkdir(exist_ok=True)

    tck_dir = ape_dir / "tcks"
    tck_dir.mkdir(exist_ok=True)

    inter_dir = roi_dir / "intersections"
    inter_dir.mkdir(exist_ok=True)

    # precompute full visual-area masks once
    area_masks = {
        area: extract_visual_area_mask(varea_map, area, roi_dir / f"{area}.nii.gz")
        for area in AREA_LABELS
    }

    color_map_list = [
        make_smooth_colormap(c.strip(), name=f"ecc_smooth_{i}")
        for i, c in enumerate(color_map.split(","))
    ]

    for idx, ecc in enumerate(ecc_bins):
        cmap = color_map_list[idx]

        for polar in polar_bins:
            use_polar = polar.lower() != "all"

            # --- combined mask (ecc only or ecc x polar) ---
            if use_polar:
                combined_mask = make_subject_patch_mask(
                    ecc_map, polar_map, ecc, polar, roi_dir / "ecc_polar"
                )
                tag = f"ecc{ecc}_polar{polar}"
            else:
                combined_mask = make_subject_patch_mask(
                    ecc_map, None, ecc, None, roi_dir / "ecc"
                )
                tag = f"ecc{ecc}"

            # --- ROI for each visual area inside this ecc/polar bin ---
            area_bin_rois = {}
            for area in AREA_LABELS:
                area_bin_rois[area] = intersect_masks(
                    combined_mask,
                    area_masks[area],
                    inter_dir / f"{area}_{tag}.nii.gz",
                )

            # --- explicit pairwise matrix ---
            M = np.zeros((n_areas, n_areas), dtype=float)

            for i, area_i in enumerate(AREA_LABELS):
                for j, area_j in enumerate(AREA_LABELS):
                        
                    if zero_diagonal and i == j:
                        M[i, j] = 0.0
                        continue
                    if j < i:
                        if symmetric:
                            M[i, j] = M[j, i]
                        continue    
                    roi1 = area_bin_rois[area_i]
                    roi2 = area_bin_rois[area_j]

                    tck_out = tck_dir / f"{area_i}_{area_j}_{tag}.tck"

                    if not tck_out.exists():
                        count = run_tckedit(
                            tract_tck,
                            roi1,
                            roi2,
                            tck_out,
                            ends_only=ends_only,
                            roi_order=roi_order,
                        )
                    else:
                        try:
                            count = int(
                                subprocess.check_output(
                                    ["tckinfo", str(tck_out), "-count"]
                                ).decode().split()[-1]
                            )
                        except Exception:
                            count = 0

                    a1 = mask_area(roi1, roi1)
                    a2 = mask_area(roi2, roi2)
                    area = (a1 + a2) / 2.0 if (a1 > 0 and a2 > 0) else max(a1, a2)
                    density = (count / area) if area > 0 else 0.0
                    M[i, j] = density
                    if symmetric:
                        M[j, i] = density
            out_csv = ape_dir / f"area_matrix_{tag}.csv"
            np.savetxt(out_csv, M, delimiter=",", fmt="%.6f")

            ecc_label = ecc.replace("_", "–")
            if use_polar:
                polar_label = polar.replace("_", "–")
                title = f"Area connectivity  ecc {ecc_label}°  polar {polar_label}°"
            else:
                title = f"Area connectivity  ecc {ecc_label}°"

            plot_area_matrix(
                M,
                title,
                ape_dir / f"area_matrix_{tag}",
                AREA_LABELS,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                log_scale=log_scale,
            )

def run_single_subject_matrix(
    tract_tck: Path,
    ecc_map: Path,
    polar_map: Path,
    varea_map: Path,
    ecc_bins: List[str],
    polar_bins: List[str],
    visual_area_a: str,
    visual_area_b: str,
    outdir: Path,
    ends_only: bool,
    roi_order: bool,
    color_map: str,
    log_scale: bool,
    vmin: Optional[float],
    vmax: Optional[float],
    fit_gaussian: bool,
    fit_truncated_gaussian_normalized: bool,
    make_dva_summary: bool,
    areas_per_bin: bool = False,
    area_matrix_method: str = "connectome",
    n_jobs: int = 1,
):
    outdir.mkdir(parents=True, exist_ok=True)

    # --- areas_per_bin mode ---
    if areas_per_bin:
        if area_matrix_method == "connectome":
            run_areas_per_bin_connectome(
                tract_tck=tract_tck,
                ecc_map=ecc_map,
                polar_map=polar_map,
                varea_map=varea_map,
                ecc_bins=ecc_bins,
                polar_bins=polar_bins,
                outdir=outdir,
                ends_only=ends_only,
                color_map=color_map,
                log_scale=log_scale,
                vmin=vmin,
                vmax=vmax,
            )
        elif area_matrix_method == "pairwise":
            run_areas_per_bin_pairwise(
                tract_tck=tract_tck,
                ecc_map=ecc_map,
                polar_map=polar_map,
                varea_map=varea_map,
                ecc_bins=ecc_bins,
                polar_bins=polar_bins,
                outdir=outdir,
                ends_only=ends_only,
                roi_order=roi_order,
                color_map=color_map,
                log_scale=log_scale,
                vmin=vmin,
                vmax=vmax,
            )
        else:
            raise ValueError(
                f"Unknown area_matrix_method '{area_matrix_method}'. "
                f"Valid options are: 'connectome', 'pairwise'."
            )
        return

    roi_dir = outdir / "ROIs"
    roi_dir.mkdir(exist_ok=True)

    # area masks
    areaA_mask = extract_visual_area_mask(varea_map, visual_area_a, roi_dir / f"{visual_area_a}.nii.gz")
    areaB_mask = extract_visual_area_mask(varea_map, visual_area_b, roi_dir / f"{visual_area_b}.nii.gz")

    M = np.zeros((len(ecc_bins), len(ecc_bins)), dtype=float)

    for i, e1 in enumerate(ecc_bins):
        for j, e2 in enumerate(ecc_bins):
            total_density = 0.0

            for polar in polar_bins:
                if polar.lower() == "all":
                    roi1_tmp = make_subject_patch_mask(ecc_map, None, e1, None, roi_dir / "ecc")
                    roi2_tmp = make_subject_patch_mask(ecc_map, None, e2, None, roi_dir / "ecc")
                    polar_tag = "all"
                else:
                    roi1_tmp = make_subject_patch_mask(ecc_map, polar_map, e1, polar, roi_dir / "ecc_polar")
                    roi2_tmp = make_subject_patch_mask(ecc_map, polar_map, e2, polar, roi_dir / "ecc_polar")
                    polar_tag = polar

                roi1 = intersect_masks(roi1_tmp, areaA_mask, roi_dir / "intersections" / f"{visual_area_a}_ecc{e1}_polar{polar_tag}.nii.gz")
                roi2 = intersect_masks(roi2_tmp, areaB_mask, roi_dir / "intersections" / f"{visual_area_b}_ecc{e2}_polar{polar_tag}.nii.gz")

                tck_out = outdir / "tcks" / f"{visual_area_a}_{visual_area_b}_{e1}_{e2}_polar{polar_tag}.tck"
                if not tck_out.exists():
                    count = run_tckedit(tract_tck, roi1, roi2, tck_out, ends_only=ends_only, roi_order=roi_order)
                else:
                    try:
                        count = int(subprocess.check_output(["tckinfo", str(tck_out), "-count"]).decode().split()[-1])
                    except Exception:
                        count = 0

                a1 = mask_area(roi1, roi1)
                a2 = mask_area(roi2, roi2)
                area = (a1 + a2) / 2.0 if (a1 > 0 and a2 > 0) else max(a1, a2)
                density = (count / area) if area > 0 else 0.0
                total_density += density

            M[i, j] = total_density

    np.savetxt(outdir / "matrix.csv", M, delimiter=",", fmt="%.6f")
    plot_matrix(
        M,
        f"Retinotopic connectivity {visual_area_a}→{visual_area_b}",
        outdir / "matrix",
        ecc_bins,
        cmap=color_map,
        vmin=vmin,
        vmax=vmax,
        log_scale=log_scale,
    )

    if make_dva_summary:
        dva_dir = outdir
        dva_dir.mkdir(exist_ok=True)
        shell_vals = compute_shell_vals(M)
        plot_dva_bar(
            shell_vals,
            out_png=dva_dir / "dva_bar.png",
            bar_x_dim=95,
            base_width=4,
            y_lim=None,
            y_decimals=4,
        )
        plot_radial_shells(M, out_png=dva_dir / "dva_radial_shells.png")

    if fit_gaussian:
        run_all_gaussian_fitting(M, ecc_bins, outdir / "gaussian_fits", fit_truncated_gaussian_normalized)
