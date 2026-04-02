# Retinotopic Connectivity Matrix (app-retinotopic-connectivity)

`app-retinotopic-connectivity` computes **single-subject retinotopic structural connectivity matrices** using a tractogram (`.tck`) and pRF-derived volumetric maps.

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

Set:

```json
"areas_per_bin": false
```

For each pair of eccentricity bins *(i, j)* the app:

1. Creates eccentricity (and optionally polar) ROI patches for two visual areas.
2. Filters streamlines using `tckedit`.
3. Computes streamline counts and normalizes by ROI volume.
4. Builds an **ecc×ecc connectivity matrix**.

---

### Mode B – `areas_per_bin` (per-bin 12×12 area matrices)

Set:

```json
"areas_per_bin": true
```

For each bin (eccentricity and optionally polar), the app computes a **12×12 visual-area connectivity matrix** across all Benson areas.

Two backends are available:

---

#### Backend 1 — `"connectome"` (default)

```json
"area_matrix_method": "connectome"
```

Pipeline:

1. Build bin mask (ecc or ecc×polar)
2. Restrict `varea.nii.gz` to the mask
3. Filter streamlines so both endpoints lie inside the bin
4. Run:

```bash
tck2connectome -symmetric -zero_diagonal -assignment_end_voxels -scale_invnodevol
```

This approach is:

* fast
* MRtrix-native
* based on endpoint assignment

---

#### Backend 2 — `"pairwise"` (strict mode)

```json
"area_matrix_method": "pairwise"
```

Pipeline:

1. Build bin mask
2. Intersect mask with each visual area
3. For every area pair `(Ai, Aj)`:

   * run:

```bash
tckedit -include ROI_i -include ROI_j [-ends_only]
```

4. Count streamlines explicitly and normalize by ROI volume

This approach is:

* stricter (explicit streamline selection)
* slower
* fully controlled (`roi_order`, `ends_only`)

---

### Important difference

| Method       | Definition of connectivity    |
| ------------ | ----------------------------- |
| `connectome` | endpoint voxel assignment     |
| `pairwise`   | explicit streamline filtering |

These methods are not equivalent and may produce different results.

---

## 2. Outputs (UPDATED)

Replace Mode B outputs:

---

### Mode B (`areas_per_bin`)

Outputs depend on backend:

#### connectome

```
output/areas_per_bin_connectome/
```

#### pairwise

```
output/areas_per_bin_pairwise/
```

Contents:

* `area_matrix_*.csv` — 12×12 matrices
* `area_matrix_*.png` — heatmaps
* `ROIs/` — bin masks and intersections
* `tcks/` — intermediate streamline subsets

---

## 3. Configuration (UPDATED)

### Mode selection

| key                  | type   | default        | meaning                        |
| -------------------- | ------ | -------------- | ------------------------------ |
| `areas_per_bin`      | bool   | false          | Enable per-bin area matrices   |
| `area_matrix_method` | string | `"connectome"` | `"connectome"` or `"pairwise"` |

---

### Pairwise-specific behavior

(only applies if `"pairwise"`)

| key         | type | default | meaning                      |
| ----------- | ---- | ------- | ---------------------------- |
| `ends_only` | bool | true    | Use endpoint-only filtering  |
| `roi_order` | bool | false   | Ordered streamline filtering |

---

### Matrix structure

Pairwise backend uses:

* `zero_diagonal = True` (diagonal entries skipped)
* matrix is intrinsically symmetric

For efficiency:

* only the upper triangle is computed
* the lower triangle is reconstructed when needed

---

## 4. Color handling (NEW)

`color_map` behavior in `areas_per_bin` mode:

| input                 | behavior                                    |
| --------------------- | ------------------------------------------- |
| `""`                  | defaults to `viridis`                       |
| `"viridis"` / `"jet"` | samples from matplotlib colormap            |
| `"red"` / `"#ff0000"` | creates a smooth custom colormap            |
| `"c1,c2,...,cN"`      | explicit colors (must match number of bins) |

---

## 5. Example configs (UPDATED)

### Pairwise strict mode

```json
{
  "areas_per_bin": true,
  "area_matrix_method": "pairwise",
  "ecc_bins": "0-2,2-4,4-6,6-8",
  "polar_bins": "all",
  "ends_only": true,
  "roi_order": false,
  "color_map": "",
  "log_scale": false
}
```

---

### Connectome mode (default)

```json
{
  "areas_per_bin": true,
  "area_matrix_method": "connectome",
  "ecc_bins": "0-2,2-4,4-6,6-8",
  "polar_bins": "all"
}
```
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

