#!/bin/bash

# Path to the key file
KEY_FILE="key-file.txt"

# Read the key file line by line
while IFS=$'\t' read -r SRR NEW_NAME; do
  # Find all files matching the SRR ID
  for FILE in ${SRR}*; do
    if [[ -e "$FILE" ]]; then
      # Extract the file extension
      EXT="${FILE#$SRR}"
      # Construct the new file name
      NEW_FILE="${NEW_NAME}${EXT}"
      # Rename the file
      mv "$FILE" "$NEW_FILE"
      echo "Renamed $FILE to $NEW_FILE"
    else
      echo "No files found for $SRR, skipping..."
    fi
  done
done < "$KEY_FILE"

echo "Renaming completed!"

