#!/usr/bin/env bash
set -euo pipefail

datatype_tags=()

function join_by { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

datatype_tags_str=""
if [ ${#datatype_tags[@]} -gt 0 ]; then
    datatype_tags_str=$(join_by , "${datatype_tags[@]}")
fi

qa_entries=()

for image in ./output/*.png; do
    if [ -f "$image" ]; then
        base=$(basename "$image" .png)

        qa_entry=$(cat <<EOF
{
  "type": "image/png",
  "name": "$base",
  "base64": "$(base64 -w 0 "$image")"
}
EOF
)
        qa_entries+=("$qa_entry")
    fi
done

if [ ${#qa_entries[@]} -eq 0 ]; then
    brainlife_json='{
        "type": "error",
        "msg": "Failed to generate output image."
    }'
else
    brainlife_json=$(printf '%s\n' "${qa_entries[@]}" | paste -sd, -)
fi

cat << EOF > product.json
{
  "datatype_tags": [],
  "brainlife": [${brainlife_json}]
}
EOF
