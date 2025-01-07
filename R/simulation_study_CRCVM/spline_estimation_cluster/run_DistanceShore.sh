#!/bin/bash

# Step 1: Prompt the user for the hyperparameters file
read -p "Enter the path to the hyperparameters file (e.g., hyperparams_set1.txt): " hyperparams_file

if [[ ! -f "$hyperparams_file" ]]; then
    echo "File not found!"
    exit 1
fi


# Step 2: Ask the user for the batch size (must be <= 50)
while true; do
    read -p "Enter the batch size (maximum 50): " batch_size
    if [[ "$batch_size" -le 50 && "$batch_size" -gt 0 ]]; then
        break
    else
        echo "Invalid input! Please enter a number between 1 and 50."
    fi
done

# Step 3: Ask the user for the number of scripts to generate
read -p "Enter the number of scripts to generate: " num_scripts

# Step 4: Update the set_up R scripts with the hyperparameters file
sed -i "s|hyperparams_file.txt|$hyperparams_file|g" set_up_rect.R
sed -i "s|hyperparams_file.txt|$hyperparams_file|g" set_up_circ.R
sed -i "s|hyperparams_file.txt|$hyperparams_file|g" set_up_fjords.R

# Step 5: Update the generate_scripts.R file with the number of scripts to generate
sed -i "s/N_SCRIPTS=200/N_SCRIPTS=$num_scripts/g" generate_scripts.R

# Step 6: Generate the R scripts by running generate_scripts.R
Rscript generate_scripts.R


# Step 7: Function to submit a batch of jobs
submit_batch() {
    local folder=$1    # Folder name (rect, circ, or fjords)
    local start=$2     # Start seed number
    local end=$3       # End seed number
    local dependency_job=$4  # Dependency job ID (optional)

    for ((k=start; k<=end; k++))
    do
        # Set job name based on seed number
        job_name="${folder}_DistanceShore$k"
        script_file="${folder}/DistanceShore$k.R"

        # Create a job submission script
        echo "#!/bin/bash" > "$job_name.sh"
        echo "#OAR -n $job_name" >> "$job_name.sh"
        echo "#OAR -l /nodes=1/core=16,walltime=03:00:00" >> "$job_name.sh"
        echo "#OAR --stdout $job_name.out" >> "$job_name.sh"
        echo "#OAR --stderr $job_name.err" >> "$job_name.sh"
        echo "#OAR --project pr-whales" >> "$job_name.sh"
        echo "source /applis/site/guix-start.sh" >> "$job_name.sh"
        echo "Rscript $script_file" >> "$job_name.sh"

        # Make the script executable
        chmod +x "$job_name.sh"

        # Submit the job, optionally with dependency
        if [[ -n $dependency_job ]]; then
            job_id=$(oarsub -S "./$job_name.sh" -a $dependency_job)
        else
            job_id=$(oarsub -S "./$job_name.sh")
        fi

        echo "Submitted $job_name with Job ID: $job_id"
    done
}

# Step 8: Loop through folders and submit jobs in batches 
folders=("rect" "circ" "fjords")
batch_size=50
start_seed=1
end_seed=200

for folder in "${folders[@]}"
do
    start_k=$start_seed

    while [[ $start_k -le $end_seed ]]
    do
        batch_end=$((start_k + batch_size - 1))

        # Ensure the last batch doesn't exceed 200
        if [[ $batch_end -gt $end_seed ]]; then
            batch_end=$end_seed
        fi

        if [[ $start_k -eq $start_seed ]]; then
            # Submit the first batch without dependency
            submit_batch $folder $start_k $batch_end
        else
            # Submit subsequent batches with dependency on the previous one
            submit_batch $folder $start_k $batch_end $previous_job_id
        fi

        # Save the last submitted job ID to use as dependency for the next batch
        previous_job_id=$(oarstat -u $(whoami) | tail -1 | awk '{print $1}')

        # Update for the next batch
        start_k=$((batch_end + 1))

        # Optionally, sleep to avoid overwhelming the scheduler
        sleep 2
    done
done

# Step 4: Run get_results.R after all jobs are completed
final_dependency_job=$(oarstat -u $(whoami) | tail -1 | awk '{print $1}')
oarsub -a $final_dependency_job -S "Rscript get_results.R"

