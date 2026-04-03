#!/usr/bin/env bash
set -euo pipefail

datatype_tags=()

function join_by { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

datatype_tags_str=""
if [ ${#datatype_tags[@]} -gt 0 ]; then
    datatype_tags_str=$(join_by , "${datatype_tags[@]}")
fi

qa_entries=()

while IFS= read -r -d '' image; do
    base="$(basename "$image" .png)"

    qa_entry="$(cat <<EOF
{
    "type": "image/png",
    "name": "$base",
    "base64": "$(base64 -w 0 "$image")"
}
EOF
)"
    qa_entries+=("$qa_entry")

done < <(find ./output -type f -name "*.png" -print0 | sort -z)

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
    "datatype_tags": [${datatype_tags_str}],
    "brainlife": [${brainlife_json}]
}
EOF
