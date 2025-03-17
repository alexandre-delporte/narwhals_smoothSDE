#!/bin/bash

echo "Files from the previous set up will be deleted. Make sure all previous calculations are over before continuing."

while true; do
    read -p "Do you want to continue? (yes/no): " answer
    case $answer in
        [Yy][Ee][Ss] ) 
            break
            ;;
        [Nn][Oo] ) 
            exit
            ;;
        * ) 
            echo "Invalid input. Please type 'yes' or 'no'."
            ;;
    esac
done

# Step 1: Ask the user for the domain of the SDE
read -p "Enter the domain for the SDE (options: rect, circ, fjords): " DOMAIN

# Validate DOMAIN input
if [[ "$DOMAIN" != "rect" && "$DOMAIN" != "circ" && "$DOMAIN" != "fjords" ]]; then
  echo "Invalid domain specified. Please use one of the following: rect, circ, fjords."
  exit 1
fi

# Verify the folder for the selected domain exists
if [[ ! -d "$DOMAIN" ]]; then
  echo "Folder for domain '$DOMAIN' not created. Check the README file for more information."
  exit 1
fi
 
# Verify the base script file exists
BASE_SCRIPT="base_script.R"

if [[ ! -f "$BASE_SCRIPT" ]]; then
  echo "File '$BASE_SCRIPT' not found."
  exit 1
fi


# Step 2: Ask the user for the number of simulations
read -p "Enter the number of simulations (e.g., 100): " NSIM

# Validate NSIM input
if ! [[ "$NSIM" =~ ^[0-9]+$ ]]; then
  echo "Invalid number of simulations. Please enter a positive integer."
  exit 1
fi

# Step 3: Modify the "generate_scripts.R" file
GENERATE_SCRIPTS="generate_scripts.R"
if [[ ! -f $GENERATE_SCRIPTS ]]; then
  echo "$GENERATE_SCRIPTS not found. Ensure the script is in the current directory."
  exit 1
fi

# Replace the value for N_SCRIPTS and domain in generate_scripts.R
sed -i "s/^N_SCRIPTS=.*/N_SCRIPTS=$NSIM/" "$GENERATE_SCRIPTS"
sed -i "s/^domain=.*/domain=\"$DOMAIN\"/" "$GENERATE_SCRIPTS"

# Step 5: Delete old R scripts, old submission bash scripts, .out and .err files and run the generate_scripts.R script locally
echo "Deleting old R scripts and log files in $DOMAIN/Rscripts..."
rm -f "$DOMAIN/Rscripts/"*.R
rm -f "$DOMAIN/Rscripts/"*.err
rm -f "$DOMAIN/Rscripts/"*.out
rm -f "$DOMAIN/Rscripts/"*.sh

echo "Generating new R scripts in $DOMAIN/Rscripts..."
Rscript "$GENERATE_SCRIPTS"


# Step 6: Ask the user for the hyperparameters file
read -p "Enter the name of the hyperparameters file (located in folder $DOMAIN): " HYPERPARAMETERS_FILE

# Validate that the hyperparameters file exists
if [[ ! -f "$DOMAIN/$HYPERPARAMETERS_FILE" ]]; then
  echo "Hyperparameters file '$HYPERPARAMETERS_FILE' not found in folder '$DOMAIN'."
  exit 1
fi

# Step 7: Modify the "set_up_DOMAIN.R" file
SETUP_SCRIPT="$DOMAIN/set_up_${DOMAIN}.R"
if [[ ! -f $SETUP_SCRIPT ]]; then
  echo "$SETUP_SCRIPT not found. Ensure the file exists in the folder '$DOMAIN'."
  exit 1
fi

# Replace the hyperparams_file line in set_up_DOMAIN.R
sed -i "s|^hyperparams_file=.*|hyperparams_file=\"$HYPERPARAMETERS_FILE\"|" "$SETUP_SCRIPT"

echo "Preparation steps completed. You can now proceed with executing the generated scripts on the cluster."

