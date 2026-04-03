# Retinotopic Connectivity Matrix (app-retinotopic-connectivity)

`app-retinotopic-connectivity` computes **single-subject retinotopic structural connectivity matrices** using a tractogram (`.tck`) and pRF-derived volumetric maps.

---

## Author

**Gabriele Amorosino**
Email: [gabriele.amorosino@utexas.edu](mailto:gabriele.amorosino@utexas.edu)

---

## Usage

### Running on Brainlife.io

Run the app from the Brainlife UI or CLI.

Notes:

* The UI may expose only a subset of advanced options
* Advanced options can be set via `config.json`

---

### Running locally

#### Prerequisites

* Singularity / Apptainer
* Access to: `docker://gamorosino/tract_align:latest`

#### Steps

```bash
git clone https://github.com/gamorosino/app-retinotopic-connectivity.git
cd app-retinotopic-connectivity
chmod +x main
./main
```

Or:

```bash
CONFIG=/path/to/config.json ./main
```

---

## Inputs

### Expected structure (`./input/`)

**Tractogram**

* `input/tracts/tractogram.tck`

**pRF maps**

* `input/prf/eccentricity.nii.gz`
* `input/prf/polarAngle.nii.gz`
* `input/prf/varea.nii.gz`

**Optional**

* `input/prf/rfWidth.nii.gz` (unused)

---

## Assumptions

1. Tractogram and maps are in the **same space**
2. `varea.nii.gz` uses Benson labels:

| label | value |
| ----- | ----: |
| V1    |     1 |
| V2    |     2 |
| V3    |     3 |
| hV4   |     4 |
| VO1   |     5 |
| VO2   |     6 |
| LO1   |     7 |
| LO2   |     8 |
| TO1   |     9 |
| TO2   |    10 |
| V3b   |    11 |
| V3a   |    12 |

---

## Execution modes

### Mode selection

```json
"mode": "area_by_area" | "bin_by_bin"
```

---

### `bin_by_bin`

Compute connectivity between eccentricity bins across two visual areas.

Pipeline:

1. Build ROI patches (eccentricity and optional polar)
2. Filter streamlines (`tckedit`)
3. Count streamlines and normalize by ROI volume
4. Build ecc × ecc matrix

Output:

* Single ecc × ecc matrix

---

### `area_by_area`

Compute connectivity across all visual areas for each bin.

Output:

* One 12 × 12 matrix per (ecc, polar)

---

### Global mode

```json
"areas_global": true
```

* Collapse bins into a single 12 × 12 matrix
* Only valid with `mode = "area_by_area"`

---

## Backends

### `connectome` (default)

Uses:

```bash
tck2connectome -symmetric -zero_diagonal -assignment_end_voxels -scale_invnodevol
```

Characteristics:

* Fast
* Endpoint-based
* MRtrix-native

---

### `pairwise`

Uses:

```bash
tckedit -include ROI_i -include ROI_j [-ends_only]
```

Characteristics:

* Explicit streamline filtering
* Slower but stricter
* Supports `ends_only` and `roi_order`

---

## Important differences

| Method     | Definition                    |
| ---------- | ----------------------------- |
| connectome | endpoint voxel assignment     |
| pairwise   | explicit streamline filtering |

These methods are not equivalent.

---

## Outputs

### `bin_by_bin`

```
output/matrix.csv
output/matrix.png
output/ROIs/
output/tcks/
```

---

### `area_by_area`

```
output/area_matrix_*.csv
output/area_matrix_*.png
output/ROIs/
output/tcks/
```

---

### Global mode

```
output/area_matrix.csv
output/area_matrix.png
```

---

## Configuration

### Mode selection

| key                | type   | default    | description                              |
| ------------------ | ------ | ---------- | ---------------------------------------- |
| mode               | string | bin_by_bin | Execution mode                           |
| areas_global       | bool   | false      | Single global matrix (area_by_area only) |
| area_matrix_method | string | connectome | connectome or pairwise                   |

---

### Visual areas (bin_by_bin only)

| key           | type   | description |
| ------------- | ------ | ----------- |
| visual_area_a | string | Source area |
| visual_area_b | string | Target area |

---

### ROI binning

| key        | type   | description            |
| ---------- | ------ | ---------------------- |
| ecc_bins   | string | "0-2,2-4,..." or "all" |
| polar_bins | string | "all" or "0-90,..."    |

---

### MRtrix filtering

| key       | type | description                          |
| --------- | ---- | ------------------------------------ |
| ends_only | bool | Endpoint-only vs full streamline     |
| roi_order | bool | Ordered endpoint matching (pairwise) |

---

### Parallelization

| key    | type | description    |
| ------ | ---- | -------------- |
| n_jobs | int  | -1 = all cores |

---

### Plot options

| key       | type   | description       |
| --------- | ------ | ----------------- |
| log_scale | bool   | Log scaling       |
| color_map | string | Colormap or color |
| range     | string | "vmin,vmax"       |

---

### Post-analysis

| key                               | type | description          |
| --------------------------------- | ---- | -------------------- |
| make_dva_summary                  | bool | Eccentricity summary |
| fit_gaussian                      | bool | Gaussian fitting     |
| fit_truncated_gaussian_normalized | bool | Normalized fit       |

---

## Color handling

### area_by_area (non-global)

* One color per eccentricity bin
* Accepts:

  * colormap
  * single color
  * list of colors

### area_by_area (global) and bin_by_bin

* Single colormap or color only

---

## Examples

### bin_by_bin

```json
{
  "mode": "bin_by_bin",
  "visual_area_a": "V1",
  "visual_area_b": "V2",
  "ecc_bins": "0-2,2-4,4-6,6-8",
  "polar_bins": "all",
  "ends_only": true,
  "roi_order": true
}
```

---

### area_by_area

```json
{
  "mode": "area_by_area",
  "area_matrix_method": "connectome",
  "ecc_bins": "0-2,2-4,4-6,6-8",
  "polar_bins": "all"
}
```

---

### global

```json
{
  "mode": "area_by_area",
  "areas_global": true,
  "area_matrix_method": "pairwise"
}
```

---

## Container execution

```bash
singularity exec -e \
  docker://gamorosino/tract_align:latest \
  micromamba run -n tract_align python3 main.py
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

