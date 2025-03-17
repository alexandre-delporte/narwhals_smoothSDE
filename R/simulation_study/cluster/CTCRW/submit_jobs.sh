#!/bin/bash

# Extract the current domain from generate_scripts.R
DOMAIN=$(grep '^domain=' generate_scripts.R | cut -d '=' -f2 | tr -d '\"')
# Extract the number of Rscripts to run
NSCRIPTS=$(grep '^N_SCRIPTS=' generate_scripts.R | cut -d '=' -f2 | tr -d '\"')
# Extract the name of the current hyperparams file
SETUP_FILE="$DOMAIN/set_up_${DOMAIN}_${TYPE}.R"
HYPERPARAMS=$(grep '^hyperparams_file=' "$SETUP_FILE" | cut -d '=' -f2 | tr -d '\"')
echo "Current domain is $DOMAIN, current simulation study type is $TYPE and current hyperparameters file is $HYPERPARAMS"

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

# Ask the user for the number of cores
read -p "Enter the number of cores for each job (max 64) " CORES
# Ask the user for the waltime
read -p "Enter the walltime (max 24:00:00 for 24h) "  WALLTIME
# Ask the user for the range of scripts to execute
read -p "There are $NSCRIPTS Rscripts to run. Enter the range of the scripts you want to run (e.g., 1-50, max 50 scripts)" RANGE

# Extract kmin and kmax from the input
kmin=$(echo "$RANGE" | cut -d '-' -f1)
kmax=$(echo "$RANGE" | cut -d '-' -f2)



# Check if DOMAIN is not empty
if [[ -z "$DOMAIN" ]]; then
    echo "Error: Unable to determine the current domain. Please set the domain in generate_scripts.R."
    exit 1
fi

# Loop through values of the seed k
for ((k=kmin; k<=kmax; k++))
do
    # Set job name
    job_name="${DOMAIN}$k"

    # Define the script path
    script_path="${DOMAIN}/Rscripts/${DOMAIN}$k.R"

    # Check if the R script exists
    if [[ ! -f "$script_path" ]]; then
        echo "Error: R script $script_path not found. Skipping job $job_name."
        continue
    fi

    # Define the folder where the job script should be created
    job_script_dir=$(dirname "$script_path")
    job_script_path="$job_script_dir/$job_name.sh"

    # Create the job submission script
    echo "#!/bin/bash" > "$job_script_path"
    echo "#OAR -n $job_name" >> "$job_script_path"
    echo "#OAR -l /nodes=1/core=$CORES,walltime=$WALLTIME" >> "$job_script_path"
    echo "#OAR --stdout $job_script_dir/$job_name.out" >> "$job_script_path"
    echo "#OAR --stderr $job_script_dir/$job_name.err" >> "$job_script_path"
    echo "#OAR --project pr-whales" >> "$job_script_path"
    echo "source /applis/site/guix-start.sh" >> "$job_script_path"
    echo "Rscript $script_path" >> "$job_script_path"

    # Make the script executable
    chmod +x "$job_script_path"

    # Submit the job
    oarsub -S "$job_script_path"

done

