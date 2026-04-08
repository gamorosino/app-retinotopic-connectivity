# Retinotopic Connectivity Matrix

`app-retinotopic-connectivity`

This app computes **retinotopic structural connectivity matrices** from a tractogram and population receptive field (pRF) maps.

The connectivity is derived from  **eccentricity**, **polar angle**, and **visual area** maps and quantified using streamline counts from an MRtrix tractogram.

---

# Author

**Gabriele Amorosino** ([gabriele.amorosino@utexas.edu](mailto:gabriele.amorosino@utexas.edu))

---

# Overview

The app measures structural connectivity between retinotopic regions of the visual cortex.

It supports two main workflows:

| Mode           | Description                                                    |
| -------------- | -------------------------------------------------------------- |
| `bin_by_bin`   | Connectivity between eccentricity bins across two visual areas |
| `area_by_area` | Connectivity between visual areas within retinotopic bins      |

Two backend implementations are available:

| Backend      | Description                                  |
| ------------ | -------------------------------------------- |
| `connectome` | Uses MRtrix `tck2connectome`                 |
| `pairwise`   | Explicit streamline filtering with `tckedit` |

---

# Inputs

Expected folder structure:

```
input/
├── tracts/
│   └── tractogram.tck
└── prf/
    ├── eccentricity.nii.gz
    ├── polarAngle.nii.gz
    └── varea.nii.gz
```

---

# Requirements

* MRtrix3
* Singularity or Apptainer
* Docker image:

```
docker://gamorosino/tract_align:latest
```

---

# Running the app

## Brainlife

Run through the Brainlife interface.

The UI exposes common options while advanced parameters can be specified via `config.json`.

---

## Running locally

Clone the repository:

```bash
git clone https://github.com/gamorosino/app-retinotopic-connectivity.git
cd app-retinotopic-connectivity
chmod +x main
./main
```

Or provide a configuration file:

```bash
CONFIG=config.json ./main
```

---

# Visual area labels

`varea.nii.gz` must use the following Benson labels:

| Label | Value |
| ----- | ----- |
| V1    | 1     |
| V2    | 2     |
| V3    | 3     |
| hV4   | 4     |
| VO1   | 5     |
| VO2   | 6     |
| LO1   | 7     |
| LO2   | 8     |
| TO1   | 9     |
| TO2   | 10    |
| V3b   | 11    |
| V3a   | 12    |

---

# Execution modes

## bin_by_bin

Computes connectivity between **eccentricity bins** across two visual areas.

Pipeline:

1. Build retinotopic ROI patches
2. Filter streamlines using MRtrix
3. Count streamlines
4. Normalize by ROI volume
5. Generate an eccentricity × eccentricity matrix

Output:

```
ecc_bins × ecc_bins matrix
```

---

## area_by_area

Computes connectivity across all visual areas for each retinotopic bin.

Each bin produces a **12 × 12 matrix** representing connectivity between visual areas.

Output:

```
one matrix per eccentricity (and optional polar) bin
```

---

## Global mode

```
"areas_global": true
```

This collapses bins into a **single visual-area connectivity matrix**.

Only valid when:

```
mode = "area_by_area"
```

---

# Connectivity backends

## connectome (default)

Uses MRtrix:

```
tck2connectome
```

Parameters:

```
-symmetric
-zero_diagonal
-assignment_end_voxels
-scale_invnodevol
```

Characteristics:

* fast
* endpoint-based assignment
* MRtrix native

---

## pairwise

Uses explicit streamline filtering:

```
tckedit -include ROI_i -include ROI_j
```

Characteristics:

* explicit streamline filtering
* slower but more controlled
* supports directional filtering

Additional options:

```
ends_only
roi_order
```

---

# Important methodological differences

| Method     | Connectivity definition       |
| ---------- | ----------------------------- |
| connectome | endpoint voxel assignment     |
| pairwise   | explicit streamline filtering |

Results from these methods **are not identical**.

---

# Configuration parameters

## Mode selection

| Parameter          | Type   | Description                    |
| ------------------ | ------ | ------------------------------ |
| mode               | string | `bin_by_bin` or `area_by_area` |
| areas_global       | bool   | compute single global matrix   |
| area_matrix_method | string | `connectome` or `pairwise`     |

---

## Visual areas (bin_by_bin only)

| Parameter     | Description        |
| ------------- | ------------------ |
| visual_area_a | source visual area |
| visual_area_b | target visual area |

---

## ROI binning

| Parameter  | Description                       |
| ---------- | --------------------------------- |
| ecc_bins   | eccentricity bins (`0_2,2_4,...`) |
| polar_bins | polar bins or `all`               |

---

## Streamline filtering

| Parameter | Description                    |
| --------- | ------------------------------ |
| ends_only | test streamline endpoints only |
| roi_order | enforce directional filtering  |

---

## Parallelization

| Parameter | Description                            |
| --------- | -------------------------------------- |
| n_jobs    | number of CPU cores (`-1` = all cores) |

---

## Plot options

| Parameter | Description             |
| --------- | ----------------------- |
| log_scale | log scale visualization |
| color_map | colormap or base color  |
| range     | color scaling limits    |

---

## Post-analysis

| Parameter                         | Description                          |
| --------------------------------- | ------------------------------------ |
| make_dva_summary                  | generate radial connectivity summary |
| fit_gaussian                      | fit Gaussian connectivity model      |
| fit_truncated_gaussian_normalized | normalized truncated Gaussian fit    |

---

# Outputs

The app exports two Brainlife-compatible datatypes.

---

# Figures datatype

```
figures/
├── images.json
└── images/
```

All generated plots are stored in `images/`.

Example:

```
figures/images/
├── matrix.png
├── dva_bar.png
├── gaussian_fit_profile.png
```

`images.json` describes each figure.

Example:

```json
{
  "images": [
    {
      "filename": "images/matrix.png",
      "name": "connectivity_matrix",
      "desc": "Retinotopic connectivity matrix"
    }
  ]
}
```

---

# Matrices datatype

```
matrices/
├── index.json
├── label.json
└── csv/
```

`csv/` contains connectivity matrices.

Example:

```
matrices/csv/
├── matrix.csv
├── area_matrix_ecc0_2.csv
```

---

### index.json

Describes matrices stored in `csv/`.

Example:

```json
[
  {
    "filename": "matrix.csv",
    "unit": "streamline_density",
    "name": "retinotopic connectivity matrix",
    "desc": "Connectivity between retinotopic bins"
  }
]
```

---

### label.json

Describes rows and columns of the matrices.

Example:

```json
[
  { "name": "self-loop", "desc": "Diagonal entries" },
  { "name": "V1", "label": "1", "voxel_value": 1 },
  { "name": "V2", "label": "2", "voxel_value": 2 }
]
```

---

# Internal working directory

Intermediate files are written to:

```
output/_work/
```

This includes:

```
ROIs/
tcks/
intersections/
```

These files are **temporary computation artifacts** and are not exported.

---

# Container execution

Example:

```bash
singularity exec -e \
  docker://gamorosino/tract_align:latest \
  micromamba run -n tract_align python3 main.py
```

---

# Data provenance

* All computations are performed in **subject space**.
* No atlas or template registration is required.
* Tractogram and pRF maps must share the **same voxel space**.

---

# Citation

If you use this app in your research, please cite:

- Hayashi, S., ... & Pestilli, F. (2024). brainlife. io: A decentralized and open-source cloud platform to support neuroscience research. Nature methods, 21(5), 809-813.
[https://doi.org/10.1038/s41592-024-02237-2](https://doi.org/10.1038/s41592-024-02237-2)

- Tournier, J. D., Smith, R., Raffelt, D., Tabbara, R., Dhollander, T., Pietsch, M., ... & Connelly, A. (2019). MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. Neuroimage, 202, 116137.
[https://doi.org/10.1016/j.neuroimage.2019.116137](https://doi.org/10.1016/j.neuroimage.2019.116137)
