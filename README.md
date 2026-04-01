# Retinotopic Connectivity Matrix (app-retinotopic-connectivity)

`app-retinotopic-connectivity` computes **single-subject retinotopic structural connectivity matrices** using an MRtrix tractogram (`.tck`) and pRF-derived volumetric maps.

Two execution modes are supported:

- **Mode A – area-to-area ecc×ecc matrix** (default): builds eccentricity/polar ROI patches for two user-specified visual areas and computes streamline density between every pair of eccentricity bins.
- **Mode B – `areas_per_ecc`**: for each eccentricity bin (and optionally each polar-angle bin), computes a **12×12 visual-area connectivity matrix** (all Benson areas) using MRtrix `tck2connectome`.

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

## Execution modes

### Mode A – area-to-area ecc×ecc matrix (default)

Set `areas_per_ecc: false` (or omit the key).

For each pair of eccentricity bins *(i, j)* the app:
1. Creates eccentricity (and optionally polar) ROI patches for each of the two visual areas.
2. Filters the tractogram with `tckedit` (`-include roi1 -include roi2 [-ends_only]`).
3. Counts the resulting streamlines and normalises by ROI volume (mm³) to obtain a **streamline density**.
4. Sums over all polar bins (if polar bins are provided).

**Outputs** (under `output/`):
- `matrix.csv` — *N×N* streamline-density matrix
- `matrix.png`

### Mode B – `areas_per_ecc` per-eccentricity 12×12 area matrix

Set `areas_per_ecc: true`.

For each eccentricity bin (and each polar bin, or all-polar if `polar_bins: "all"`) the app:
1. Creates a combined mask (eccentricity only, or eccentricity ∩ polar).
2. Multiplies `varea.nii.gz` by the combined mask to obtain a **label image** restricted to that bin.
3. Filters the tractogram with `tckedit` so that **both streamline endpoints** lie within the combined mask (`-include mask -include mask [-ends_only]`).
4. Runs `tck2connectome` with flags `-symmetric -zero_diagonal -assignment_end_voxels -scale_invnodevol -force` to compute the 12×12 area connectivity matrix.
5. Saves the result as a CSV and plots a PNG.

**Output naming**:
- Eccentricity-only (`polar_bins: "all"`): `area_matrix_ecc{bin}.csv` / `.png`
- With polar bins: `area_matrix_ecc{bin}_polar{bin}.csv` / `.png`

**Outputs** (under `output/areas_per_ecc/`):
- `area_matrix_ecc*.csv` — 12×12 connectivity matrices
- `area_matrix_ecc*.png` — heatmap plots
- `ROIs/` — generated eccentricity/polar/label masks
- `tcks/` — filtered `.tck` files

---

## Outputs

### Mode A
- `output/matrix.csv`
- `output/matrix.png`
- `output/ROIs/` — generated masks
- `output/tcks/` — filtered `.tck` files
- `output/gaussian_fits/` *(optional)*
- `output/dva_summary/` *(optional)*

### Mode B
- `output/areas_per_ecc/area_matrix_ecc*.csv`
- `output/areas_per_ecc/area_matrix_ecc*.png`
- `output/areas_per_ecc/ROIs/`
- `output/areas_per_ecc/tcks/`

---

## Configuration (`config.json`)

All keys are **snake_case**. The wrapper reads `config.json` and passes options to Python.

### 1) Mode selection

| config key | type | default | meaning |
|---|---|---|---|
| `areas_per_ecc` | bool | `false` | Run Mode B (12×12 per-ecc area matrix) |

### 2) Visual areas (Mode A only)

| config key | type | meaning |
|---|---|---|
| `visual_area_a` | string | Source area (e.g. `"V1"`) |
| `visual_area_b` | string | Target area (e.g. `"V2"`) |

### 3) ROI binning

| config key | type | meaning |
|---|---|---|
| `ecc_bins` | string | Comma-separated bins: `"0-2,2-4,4-6,6-8,8-90"` |
| `polar_bins` | string | `"all"` or comma-separated bins `"0-15,30-60,..."` |

Notes:
- You can use either `-` or `_` in bin strings; the app normalizes them internally.
- If `polar_bins` is `"all"`, no polar restriction is applied (eccentricity-only ROIs).

### 4) MRtrix filtering options

| config key | type | meaning |
|---|---|---|
| `ends_only` | bool | Use `tckedit -ends_only` |
| `roi_order` | bool | Ordered ROI pairs (Mode A only) |

### 5) Plot options

| config key | type | meaning |
|---|---|---|
| `log_scale` | bool | Log color scaling for the matrix heatmap |
| `color_map` | string | Matplotlib colormap |
| `range` | string | Optional `"vmin,vmax"` (e.g. `"0,0.02"`). Empty = auto |

### 6) Optional post-analysis (Mode A only)

| config key | type | meaning |
|---|---|---|
| `fit_gaussian` | bool | Fit Gaussian models to each row of the matrix |
| `fit_truncated_gaussian_normalized` | bool | Normalize truncated Gaussian |
| `make_dva_summary` | bool | Create DVA shell summary plots |

---

## Example `config.json` – Mode A (default)

```json
{
  "areas_per_ecc": false,
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

## Example `config.json` – Mode B (`areas_per_ecc`)

```json
{
  "areas_per_ecc": true,
  "ecc_bins": "0-2,2-4,4-6,6-8",
  "polar_bins": "all",
  "ends_only": true,
  "log_scale": false,
  "color_map": "viridis",
  "range": ""
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
- **MRtrix3**: J.-D. Tournier, et al. (2019). MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. *NeuroImage*, 202. DOI: 10.1016/j.neuroimage.2019.116137

---

## Data provenance notes

- All computation is performed in subject space using the provided pRF maps.
- No atlas or template is required; the app uses only the subject-specific `eccentricity.nii.gz`, `polarAngle.nii.gz`, and `varea.nii.gz`.
- Tractogram and pRF maps must be in the same voxel space (same affine).

