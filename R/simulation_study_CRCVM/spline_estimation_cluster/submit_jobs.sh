#!/bin/bash

# Extract the current domain from generate_scripts.R
DOMAIN=$(grep '^domain=' generate_scripts.R | cut -d '=' -f2 | tr -d '\"')

# Check if DOMAIN is not empty
if [[ -z "$DOMAIN" ]]; then
    echo "Error: Unable to determine the current domain. Please set the domain in generate_scripts.R."
    exit 1
fi

# Loop through values of k from 1 to 50
for ((k=1; k<=50; k++))
do
    # Set job name
    job_name="${DOMAIN}_Ishore$k"

    # Navigate to the correct domain subfolder
    script_path="${DOMAIN}/Rscripts/${DOMAIN}_DistanceShore$k.R"

    # Check if the R script exists
    if [[ ! -f "$script_path" ]]; then
        echo "Error: R script $script_path not found. Skipping job $job_name."
        continue
    fi

    # Create the job submission script
    echo "#!/bin/bash" > "$job_name.sh"
    echo "#OAR -n $job_name" >> "$job_name.sh"
    echo "#OAR -l /nodes=1/core=16,walltime=03:00:00" >> "$job_name.sh"
    echo "#OAR --stdout $job_name.out" >> "$job_name.sh"
    echo "#OAR --stderr $job_name.err" >> "$job_name.sh"
    echo "#OAR --project pr-whales" >> "$job_name.sh"
    echo "source /applis/site/guix-start.sh" >> "$job_name.sh"
    echo "Rscript $script_path" >> "$job_name.sh"

    # Make the script executable
    chmod +x "$job_name.sh"

    # Submit the job
    oarsub -S "./$job_name.sh"

done

