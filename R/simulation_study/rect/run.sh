#!/bin/bash

# Loop through values of k from 1 to 100
for ((k=51; k<=100; k++))
do
    # Set job name
    job_name="rect_DistanceShore$k"

    # Create submission script  with OAR directives
    echo "#!/bin/bash" > "$job_name.sh"
    echo "#OAR -n $job_name" >> "$job_name.sh"
    echo "#OAR -l /core=16,walltime=01:30:00" >> "$job_name.sh"
 
    echo "#OAR --stdout $job_name.out" >> "$job_name.sh"
    echo "#OAR --stderr $job_name.err" >> "$job_name.sh"
    echo "#OAR --project pr-formation-ced-calcul" >> "$job_name.sh"
    echo "source /applis/site/guix-start.sh" >> "$job_name.sh"
    echo "Rscript rect_DistanceShore$k.R" >> "$job_name.sh"
    
    #Make the script executable
    chmod +x "$job_name.sh"

    # Submit the job
    oarsub -S ./"$job_name.sh"    
    
    # Optionally, you can add a delay between job submissions
    #sleep 600
done

