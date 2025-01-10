#!/bin/bash

# Retrieve the current domain
DOMAIN=$(grep '^domain=' generate_scripts.R | cut -d '=' -f2 | tr -d '\"')

# Retrieve the current hyperparameters file
HYPERPARAMETERS_FILE=$(grep '^hyperparams_file <- ' "$DOMAIN/set_up_$DOMAIN.R" | awk -F'"' '{print $2}')

# Output the results
echo "Current Domain: $DOMAIN"
echo "Current Hyperparameters File: $HYPERPARAMETERS_FILE"


