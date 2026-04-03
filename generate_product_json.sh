#!/usr/bin/env bash
set -euo pipefail

datatype_tags=()
echo "writing out ${2}"

function join_by { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

datatype_tags_str=""
if [ ${#datatype_tags[@]} -gt 0 ]; then
    datatype_tags_str=$(join_by , "${datatype_tags[@]}")
fi

mkdir -p qa/previews

qa_entries=()

while IFS= read -r -d '' image; do
    #base="$(basename "$image" .png)"
    #jpg="qa/previews/${base}.jpg"

    #convert "$image" -resize 50% -trim -quality 90 "$jpg"

    qa_entry="{
        \"type\": \"image/png\",
        \"name\": \"$base\",
        \"base64\": \"$(base64 -w 0 "$image")\"
    }"

    qa_entries+=("$qa_entry")
done < <(find "$1" -type f -name "*.png" -print0 | sort -z)

if [ ${#qa_entries[@]} -eq 0 ]; then
    qa='{ 
        "type": "error",
        "msg": "Failed to generate output image."
    }'
    brainlife_json="$qa"
else
    brainlife_json=$(printf "%s\n" "${qa_entries[@]}" | paste -sd, -)
fi

cat << EOF > "$2"
{
    "datatype_tags": [${datatype_tags_str}],
    "brainlife": [${brainlife_json}]
}
EOF
