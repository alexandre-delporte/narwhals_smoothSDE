#!/bin/bash

hyperparams_files=("hyperparams_set1.txt" "hyperparams_set2.txt" "hyperparams_set3.txt") 

for hp_file in "${hyperparams_files[@]}"; do
    job_name="job_$(basename "$hp_file" .txt).sh"
    
    # Create a unique job script for each hyperparameter set
    cat <<EOF > "$job_name"
#!/bin/bash
#OAR -n $job_name
#OAR -l /nodes=1/core=8,walltime=03:00:00
#OAR --stdout $job_name.out
#OAR --stderr $job_name.err
#OAR --project pr-whales

module load R
Rscript fit_model.R "$hp_file"
EOF
    # Make the script executable
    chmod +x "$job_name"
    # Submit the job
    oarsub -S "$job_name"
    sleep 1  
done

