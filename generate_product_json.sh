#!/usr/bin/env bash
set -euo pipefail

outdir=${1}
product_json=${2}

datatype_tags_str='"neuro/tractography","retinotopy"'

brainlife_entries=()

mapfile -t pngs < <(find "$outdir" -type f -name '*.png' | sort)

if [ ${#pngs[@]} -eq 0 ]; then
    brainlife_entries+=('{
        "type": "error",
        "msg": "Failed to generate output image."
    }')
else
    for image in "${pngs[@]}"; do
        name="$(basename "$image" .png)"

        qa="$(cat <<EOF
{
    "type": "image/png",
    "name": "$name",
    "base64": "$(base64 -w 0 "$image")"
}
EOF
)"
        brainlife_entries+=("$qa")
    done
fi

brainlife_json="$(printf '%s\n' "${brainlife_entries[@]}" | paste -sd, -)"

cat <<EOF > "$product_json"
{
    "datatype_tags": [$datatype_tags_str],
    "brainlife": [$brainlife_json]
}
EOF
