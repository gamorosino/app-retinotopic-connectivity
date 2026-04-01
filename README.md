# Retinotopic Connectivity Matrix (app-retinotopic-connectivity)

`app-retinotopic-connectivity` computes a **single-subject retinotopic structural connectivity matrix** between two visual areas (e.g., **V1 → V2**) using:
- an MRtrix tractogram (`.tck`)
- pRF-derived volumetric maps (`eccentricity.nii.gz`, `polarAngle.nii.gz`, `varea.nii.gz`)

The app builds eccentricity and (optionally) polar-angle ROI patches, intersects them with the requested visual areas, filters streamlines with `tckedit`, and computes a **streamline density** per ROI-pair normalized by ROI volume (mm³). It also produces a heatmap and can optionally run Gaussian fitting / DVA summary plots.

Repository contents:
- `main`: bash wrapper that reads `config.json` and executes the app inside `docker://gamorosino/tract_align:latest`.
- `main.py`: Python CLI entrypoint (Brainlife-friendly).
- `retinotopic_connectivity/connectivity.py`: core implementation (ROI generation, streamline filtering, matrix computation, plotting).
- `utils.py`, `fit_gaussian_connectivity.py`, `diagonal_to_braplot.py`: vendored helper modules from the VISCONTI codebase (used for binning utilities, Gaussian fitting, and DVA summary plots).

---

## Author

**Gabriele Amorosino**  
Email: gabriele.amorosino@utexas.edu

---

## Usage

### Running on Brainlife.io
Run the app from the Brainlife UI / CLI as usual.

Notes:
- Brainlife UI may expose only a subset of advanced options.
- Advanced options can be set via `config.json` when running locally (or if your platform supports custom config injection).

### Running locally

#### Prerequisites
- Singularity / Apptainer
- A system that can pull `docker://gamorosino/tract_align:latest`

#### Steps
```bash
git clone https://github.com/gamorosino/app-retinotopic-connectivity.git
cd app-retinotopic-connectivity
chmod +x main
./main
```

By default `main` reads `config.json` in the current directory.
You can also run with:
```bash
CONFIG=/path/to/config.json ./main
```

---

## Inputs

### Expected mounted paths

This app expects inputs to be available under `./input/` (Brainlife-style):

**Tractogram**
- `input/tracts/tractogram.tck`  
  (fallback: the first `input/tracts/*.tck` if `tractogram.tck` is not present)

**pRF maps** (required)
- `input/prf/eccentricity.nii.gz`
- `input/prf/polarAngle.nii.gz`
- `input/prf/varea.nii.gz`

**Optional**
- `input/prf/rfWidth.nii.gz` (currently unused by the matrix computation)

### Core assumptions (important)
1. The tractogram and the pRF maps are in the **same space** (voxel grid / affine alignment).
2. `varea.nii.gz` uses **Benson-style integer labels** with this ordering:

| label | value |
|---|---:|
| V1 | 1 |
| V2 | 2 |
| V3 | 3 |
| hV4 | 4 |
| VO1 | 5 |
| VO2 | 6 |
| LO1 | 7 |
| LO2 | 8 |
| TO1 | 9 |
| TO2 | 10 |
| V3b | 11 |
| V3a | 12 |

---

## Outputs

Outputs are written under `./output/`:

**Always produced**
- `output/matrix.csv` — (N×N) connectivity matrix (streamline density)
- `output/matrix.png` (and possibly `.pdf` / `.svg` depending on `save_figure()`)

**Intermediate files**
- `output/ROIs/` — generated eccentricity/polar/area/intersection masks
- `output/tcks/` — filtered `.tck` files for each ROI pair (can be large)

**Optional (if enabled in config)**
- `output/gaussian_fits/` — Gaussian and truncated-Gaussian row fits + plots
- `output/dva_summary/` — shell barplot + radial shells plot

---

## Configuration (`config.json`)

All keys are **snake_case**. The wrapper reads `config.json` and passes options to Python.

### 1) Visual areas

| config key | type | meaning |
|---|---|---|
| `visual_area_a` | string | Source area (e.g. `"V1"`) |
| `visual_area_b` | string | Target area (e.g. `"V2"`) |

### 2) ROI binning

| config key | type | meaning |
|---|---|---|
| `ecc_bins` | string | Comma-separated bins: `"0-2,2-4,4-6,6-8,8-90"` |
| `polar_bins` | string | `"all"` or comma-separated bins `"0-15,30-60,..."` |

Notes:
- You can use either `-` or `_` in bin strings; the app normalizes them internally.
- If `polar_bins` is `"all"`, no polar restriction is applied (eccentricity-only ROIs).

### 3) MRtrix filtering options

| config key | type | meaning |
|---|---|---|
| `ends_only` | bool | Use `tckedit -ends_only` |
| `roi_order` | bool | Ordered ROI pairs (ROI1→ROI2 and ROI2→ROI1 merged) |

### 4) Plot options

| config key | type | meaning |
|---|---|---|
| `log_scale` | bool | Log color scaling for the matrix heatmap |
| `color_map` | string | Matplotlib colormap |
| `range` | string | Optional `"vmin,vmax"` (e.g. `"0,0.02"`). Empty = auto |

### 5) Optional post-analysis

| config key | type | meaning |
|---|---|---|
| `fit_gaussian` | bool | Fit Gaussian models to each row of the matrix |
| `fit_truncated_gaussian_normalized` | bool | Normalize truncated Gaussian (only if `fit_gaussian=true`) |
| `make_dva_summary` | bool | Create DVA shell summary plots |

---

## Example `config.json`

```json
{
  "visual_area_a": "V1",
  "visual_area_b": "V2",
  "ecc_bins": "0-2,2-4,4-6,6-8,8-90",
  "polar_bins": "all",
  "ends_only": true,
  "roi_order": false,
  "log_scale": false,
  "color_map": "viridis",
  "range": "",
  "fit_gaussian": false,
  "fit_truncated_gaussian_normalized": false,
  "make_dva_summary": false
}
```

---

## Container execution

This app runs:

```bash
singularity exec -e \
  docker://gamorosino/tract_align:latest \
  micromamba run -n tract_align python3 main.py ...
```

---

## Citation

If you use this app in your research, please cite:

- **Brainlife.io**: Hayashi, S., et al. (2024). *Nature Methods, 21*(5), 809–813. DOI: 10.1038/s41592-024-02237-2
- **MRtrix3**: J.-D. Tournier, et al. (2019) MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. NeuroImage, 202
