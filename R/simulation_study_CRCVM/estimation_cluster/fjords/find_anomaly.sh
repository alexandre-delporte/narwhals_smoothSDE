#!/bin/bash

# Prefix to search for
prefix='"nu.s(ID)"'

# Iterate over the matching CSV files
for file in result_fjords_60h_12ID_1km_DistanceShore*.csv; do
    echo "Processing file: $file"
    # Use grep to find lines starting with the prefix
    grep "^$prefix" "$file"
    echo ""
done
