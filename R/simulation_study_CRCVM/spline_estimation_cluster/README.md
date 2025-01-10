Execute set_up_simulation_study.sh to choose the domain, the number of simulations and the hyperparameters file for the simulation study on the clusters.
Options for the domain are "rect","circ" and "fjords". Folders "rect","circ" and "fjords" need to be in the current folder and must contain a set_up file named "set_up_domain.R" (where domain is the name of the domain) that creates the domain geometry, read the hyperparameter file and define the smooth function for parameter omega in the SDE. Hyperparameter files are .txt files stored in folder "domain" in which the parameter for the simulation study are defined. These parameters are :
N_ID_LOW,N_ID_HIGH,TMAX,DELTA,BY_LF,BY_HF,SP_DF,TAU_0,NU_0,SIGMA_TAU,SIGMA_NU,SIGMA_OBS_HIGH,SIGMA_OBS_LOW,DMIN,DMAX,PAR0,
SIGMA_OBS0,SIGMA_TAU_0,SIGMA_NU_0,A,D0,D1,SIGMA_THETA,SIGMA_D,B,D_LOW,D_UP.
Scripts for the simulation study with the current set_up are generated in folder "/domain/Rscripts".

The jobs can then be submitted with the bashscript "submit_jobs.sh". Before submitting the jobs, the current set_up may be checked by executing "check_current_settings.sh". Maximum 50 jobs can be submitted at a time. Project name, number of cores and nodes as well as walltime can be set by changing the "submit_jobs.sh" script.
Once the calculation is done, the results of the parameter estimates are stored in .csv files in a folder "results_hyperparams" where hyperparams is the name of the hyperparameter file that was used for the simulation study.
