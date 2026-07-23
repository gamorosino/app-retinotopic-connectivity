#!/usr/bin/env bash
set -euo pipefail

# Build config.json from --<key> <value> flags, then run ./main with it.
#
# Every "--foo bar" pair is written into config.json as "foo": bar, with the
# value type auto-detected (true/false -> boolean, numeric -> number, else
# string). This mirrors the keys read by ./main and main.py, e.g.:
#
#   ./main_cli.sh \
#     --track /input/tracts/track.tck \
#     --eccentricity /input/prf/eccentricity.nii.gz \
#     --polarAngle /input/prf/polarAngle.nii.gz \
#     --varea /input/prf/varea.nii.gz \
#     --mode area_by_area \
#     --ecc_bins "0-2,2-4,4-6,6-8,8-90" \
#     --ends_only true
#
# Pass --dry-run to only generate and print config.json without executing.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_OUT="config.json"
CONFIG_OUT_SET=0
OUTDIR=""
DRY_RUN=0
NO_CONTAINER=0
declare -a JSON_FRAGMENTS=()

usage() {
  cat <<'EOF'
Usage: ./main_cli.sh --<key> <value> [--<key> <value> ...] [--config PATH] [--dry-run]

Any --<key> <value> pair becomes a "<key>": <value> entry in config.json.
A bare --<key> (no value, or followed by another --flag) is treated as
--<key> true, e.g. `--ends_only --roi_order` sets both to true.
By default, ./main is then executed with the generated config.

Common keys (see README.md for the full list):
  --track                     tractogram (.tck) [required]
  --parc                      parcellation image (enables parcellation mode)
  --eccentricity               eccentricity map (.nii.gz)
  --polarAngle                 polar angle map (.nii.gz)
  --varea                      visual area map (.nii.gz)
  --mode                        bin_by_bin | area_by_area | parcellation
  --visual_area_a / --visual_area_b
  --ecc_bins                    e.g. "0-2,2-4,4-6,6-8,8-90"
  --polar_bins                  e.g. "all"
  --ends_only / --roi_order
  --matrix_elements             eccentricity | polar
  --area_matrix_method          connectome | pairwise
  --log_scale / --color_map / --range
  --fit_gaussian / --fit_truncated_gaussian_normalized
  --make_dva_summary
  --n_jobs

Options:
  --config PATH     write generated config to PATH instead of ./config.json
                    (default becomes <OUTDIR>/config.json when --output_dir is set)
  --output_dir DIR  CLI-only: run ./main from inside DIR, so output/, figures/,
                    matrices/ and product.json land there instead of the repo
                    root.
  --workdir DIR     relocate the heavy intermediate computation directory
                    (ROIs/tcks/etc.) to DIR instead of nesting it under
                    output/_work. Useful to point scratch I/O at fast local
                    storage while --output_dir keeps final results elsewhere.
  --dry-run         only generate and print config.json; do not run ./main
  --no-container    skip Singularity entirely and run against
                    locally-installed software instead (a native MRtrix3 +
                    Python env with this app's dependencies already on
                    PATH/importable). Uses a local
                    `micromamba run -n tract_align` if that environment
                    exists, else runs commands directly; --compress falls
                    back to a local `convert`/`magick` (ImageMagick).
  -h, --help        show this help
EOF
}

to_json_fragment() {
  local key="$1" value="$2"
  # Resolve file/dir path values to absolute paths so they keep resolving
  # correctly even if we cd into --output_dir before running ./main.
  if [[ -e "$value" ]]; then
    value="$(realpath "$value")"
  fi
  jq -n --arg k "$key" --arg v "$value" '
    if ($v == "true") then {($k): true}
    elif ($v == "false") then {($k): false}
    elif ($v | test("^-?[0-9]+$")) then {($k): ($v | tonumber)}
    elif ($v | test("^-?[0-9]*\\.[0-9]+$")) then {($k): ($v | tonumber)}
    else {($k): $v}
    end
  '
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --config)
      [[ $# -ge 2 ]] || { echo "ERROR: missing value for --config" >&2; exit 1; }
      CONFIG_OUT="$2"
      CONFIG_OUT_SET=1
      shift 2
      ;;
    --output_dir)
      [[ $# -ge 2 ]] || { echo "ERROR: missing value for --output_dir" >&2; exit 1; }
      OUTDIR="$2"
      shift 2
      ;;
    --workdir)
      [[ $# -ge 2 ]] || { echo "ERROR: missing value for --workdir" >&2; exit 1; }
      mkdir -p "$2"
      JSON_FRAGMENTS+=("$(to_json_fragment "workdir" "$(realpath "$2")")")
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --no-container)
      NO_CONTAINER=1
      shift
      ;;
    --*=*)
      key="${1#--}"
      value="${key#*=}"
      key="${key%%=*}"
      JSON_FRAGMENTS+=("$(to_json_fragment "$key" "$value")")
      shift
      ;;
    --*)
      key="${1#--}"
      # Bare boolean flag (no value, or immediately followed by another
      # --flag) defaults to true, e.g. `--ends_only --roi_order`.
      if [[ $# -ge 2 && "$2" != --* ]]; then
        value="$2"
        shift 2
      else
        value="true"
        shift
      fi
      JSON_FRAGMENTS+=("$(to_json_fragment "$key" "$value")")
      ;;
    *)
      echo "ERROR: unrecognized argument '$1'" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ ${#JSON_FRAGMENTS[@]} -eq 0 ]]; then
  echo "ERROR: no --<key> <value> config options provided" >&2
  usage
  exit 1
fi

if [[ -n "${OUTDIR}" ]]; then
  mkdir -p "${OUTDIR}"
  OUTDIR="$(cd "${OUTDIR}" && pwd)"
  if [[ "${CONFIG_OUT_SET}" -eq 0 ]]; then
    CONFIG_OUT="${OUTDIR}/config.json"
  fi
fi

printf '%s\n' "${JSON_FRAGMENTS[@]}" | jq -s 'add' > "${CONFIG_OUT}"
CONFIG_OUT="$(realpath "${CONFIG_OUT}")"

echo "Wrote ${CONFIG_OUT}:"
cat "${CONFIG_OUT}"

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo
  echo "Dry run: not executing ./main."
  exit 0
fi

if [[ -n "${OUTDIR}" ]]; then
  echo
  echo "Running ./main in ${OUTDIR} (output/, figures/, matrices/, product.json will land there) ..."
  ( cd "${OUTDIR}" && CONFIG="${CONFIG_OUT}" NO_CONTAINER="${NO_CONTAINER}" bash "${SCRIPT_DIR}/main" )
else
  echo
  echo "Running ./main with CONFIG=${CONFIG_OUT} ..."
  CONFIG="${CONFIG_OUT}" NO_CONTAINER="${NO_CONTAINER}" bash "${SCRIPT_DIR}/main"
fi
