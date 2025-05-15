#!/bin/bash

for file in *.txt; do
    # Count the number of lines in the last column that are not 0
    count=$(awk '$NF != 0 {count++} END {print count}' "$file")
    # Print the filename and the count
    echo "$file: $count"
done