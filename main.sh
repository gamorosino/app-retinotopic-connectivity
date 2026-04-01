#!/bin/bash

# Check if config.json exists
CONFIG_FILE="config.json"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Configuration file not found!"
    exit 1
fi

# Parse file paths and configuration using jq
INPUT_FILE=$(jq -r '.input_file' "$CONFIG_FILE")
OUTPUT_DIR=$(jq -r '.output_directory' "$CONFIG_FILE")

# Example usage of the parsed values
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"

# Add any additional processing using the parsed values here...