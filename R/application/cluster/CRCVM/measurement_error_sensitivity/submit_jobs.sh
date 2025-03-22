#!/bin/bash

hyperparams_files=("hyperparams_set1.txt" "hyperparams_set2.txt" "hyperparams_set3.txt" "hyperparams_set4.txt" "hyperparams_set5.txt") 


# Ask the user for the number of cores
read -p "Enter the number of cores for each job (max 64) " CORES
# Ask the user for the waltime
read -p "Enter the walltime (max 24:00:00 for 24h) "  WALLTIME

for hp_file in "${hyperparams_files[@]}"; do
    job_name="job_$(basename "$hp_file" .txt).sh"
    
    # Create a unique job script for each hyperparameter set
	echo "#!/bin/bash" > "$job_name"
	echo "#OAR -n job_$(basename "$hp_file" .txt)" >> "$job_name"
	echo "#OAR -l /nodes=1/core=$CORES,walltime=$WALLTIME" >> "$job_name"
	echo "#OAR --stdout job_$(basename "$hp_file" .txt).out" >> "$job_name"
	echo "#OAR --stderr job_$(basename "$hp_file" .txt).err" >> "$job_name"
	echo "#OAR --project pr-whales" >> "$job_name"

	echo "source /applis/site/guix-start.sh" >> "$job_name"
	echo "Rscript fit.R '$hp_file'" >> "$job_name"

    	# Make the script executable
    	chmod +x "$job_name"
    	# Submit the job
    	oarsub -S "./$job_name"
    	sleep 1  
done

