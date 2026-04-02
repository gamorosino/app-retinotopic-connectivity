#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

from retinotopic_connectivity.connectivity import (
    run_single_subject_matrix,
    make_smooth_colormap,
)
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgb


def resolve_area_bin_colors(color_map: str, n_bins: int) -> str:
    color_map = (color_map or "").strip()

    # explicit list of colors
    if "," in color_map:
        colors = [c.strip() for c in color_map.split(",") if c.strip()]
        if len(colors) != n_bins:
            raise ValueError(
                f"Expected {n_bins} colors for areas_per_bin, got {len(colors)} "
                f"in color_map='{color_map}'."
            )
        for c in colors:
            try:
                to_rgb(c)
            except ValueError:
                raise ValueError(f"Invalid color '{c}' in color_map list.")
        return ",".join(colors)

    # default colormap
    if color_map == "":
        cmap = plt.get_cmap("viridis")
    else:
        # first try as matplotlib colormap name
        try:
            cmap = plt.get_cmap(color_map)
        except ValueError:
            # otherwise try as a single color
            try:
                to_rgb(color_map)
                cmap = make_smooth_colormap(color_map, name=f"custom_{color_map}")
            except ValueError:
                raise ValueError(
                    f"'{color_map}' is neither a valid matplotlib colormap name "
                    f"nor a valid color."
                )

    colors = [to_hex(cmap(i / max(n_bins - 1, 1))) for i in range(n_bins)]
    return ",".join(colors)
    
def _get(cfg, key, default=None):
    if cfg is None:
        return default
    return cfg.get(key, default)


def main():
    p = argparse.ArgumentParser(
        description="Compute single-subject retinotopic connectivity matrix (Brainlife app)."
    )
    p.add_argument("--tck", required=True, type=str, help="Input tractogram (.tck)")
    p.add_argument("--ecc", required=True, type=str, help="eccentricity.nii.gz")
    p.add_argument("--polar", required=True, type=str, help="polarAngle.nii.gz")
    p.add_argument("--varea", required=True, type=str, help="varea.nii.gz")
    p.add_argument("--outdir", required=True, type=str, help="Output directory")
    p.add_argument("--config", required=False, type=str, default=None, help="Optional config.json")
    args = p.parse_args()

    cfg = None
    if args.config:
        cfg = json.loads(Path(args.config).read_text())

    visual_area_a = _get(cfg, "visual_area_a", "V1")
    visual_area_b = _get(cfg, "visual_area_b", "V2")

    ecc_bins_str = _get(cfg, "ecc_bins", "0-2,2-4,4-6,6-8,8-90")
    polar_bins_str = _get(cfg, "polar_bins", "all")

    # normalize to underscore format used internally
    ecc_bins = [b.strip().replace("-", "_") for b in ecc_bins_str.split(",") if b.strip()]
    polar_bins = [b.strip().replace("-", "_") for b in polar_bins_str.split(",") if b.strip()]
    if len(polar_bins) == 0:
        polar_bins = ["all"]

    ends_only = bool(_get(cfg, "ends_only", True))
    roi_order = bool(_get(cfg, "roi_order", False))

    log_scale = bool(_get(cfg, "log_scale", False))
    color_map = _get(cfg, "color_map", "viridis")

    range_str = _get(cfg, "range", "")
    vmin = vmax = None
    if isinstance(range_str, str) and range_str.strip():
        parts = [x.strip() for x in range_str.split(",")]
        if len(parts) == 2:
            vmin = float(parts[0])
            vmax = float(parts[1])

    fit_gaussian = bool(_get(cfg, "fit_gaussian", False))
    fit_trunc_norm = bool(_get(cfg, "fit_truncated_gaussian_normalized", False))
    make_dva_summary = bool(_get(cfg, "make_dva_summary", False))

    # renamed from areas_per_ecc -> areas_per_bin
    areas_per_bin = bool(_get(cfg, "areas_per_bin", False))

    # new: choose backend for area-per-bin mode
    area_matrix_method = str(_get(cfg, "area_matrix_method", "connectome")).strip().lower()
    if area_matrix_method not in {"connectome", "pairwise"}:
        raise ValueError(
            f"Invalid area_matrix_method '{area_matrix_method}'. "
            f"Valid options are: 'connectome', 'pairwise'."
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if areas_per_bin:
        color_map = resolve_area_bin_colors(color_map, len(ecc_bins))
    else:
        if color_map == "":
            color_map = "hot"

    run_single_subject_matrix(
        tract_tck=Path(args.tck),
        ecc_map=Path(args.ecc),
        polar_map=Path(args.polar),
        varea_map=Path(args.varea),
        ecc_bins=ecc_bins,
        polar_bins=polar_bins,
        visual_area_a=visual_area_a,
        visual_area_b=visual_area_b,
        outdir=outdir,
        ends_only=ends_only,
        roi_order=roi_order,
        color_map=color_map,
        log_scale=log_scale,
        vmin=vmin,
        vmax=vmax,
        fit_gaussian=fit_gaussian,
        fit_truncated_gaussian_normalized=fit_trunc_norm,
        make_dva_summary=make_dva_summary,
        areas_per_bin=areas_per_bin,
        area_matrix_method=area_matrix_method,
    )


if __name__ == "__main__":
    main()
