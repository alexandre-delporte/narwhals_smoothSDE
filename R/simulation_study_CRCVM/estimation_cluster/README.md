Execute set_up_simulation_study.sh to choose the domain, the type of simulation study, the number of simulations and the hyperparameters file for the simulation study on the clusters.
Options for the domain are "rect","circ" and "fjords". Folders "rect","circ" and "fjords" need to be in the current folder.
Current options for the type of simulation study are :
- "Ishore", in which the parameter omega is estimated using tensor splines of the covariates thete and Ishore, the parameters tau and nu are constant with random effects and the measurement error is fixed. Covariates are computed from the noisy observations as in real life applications.
- "parametric", in which the smooth function omega is fixed to its true value, the parameters tau and nu are constant WITHOUT random effects and the measurement error is fixed. Covariates are computed from true position and velocity.
Folders "rect", "fjords" and "circ" must contain a set_up file named "set_up_domain_type.R" (where domain is the name of the domain, and type is the type of simulation study) that creates the domain geometry, read the hyperparameter file and define the smooth function for parameter omega in the SDE.
Hyperparameter files are .txt files stored in folder "domain" in which the parameter for the simulation study are defined. These parameters might differ from one type of simulation study to another. For the type "Ishore", they are :
N_ID_LOW,N_ID_HIGH,TMAX,DELTA,BY_LF,BY_HF,SP_DF,TAU_0,NU_0,SIGMA_TAU,SIGMA_NU,SIGMA_OBS_HIGH,SIGMA_OBS_LOW,DMIN,DMAX,PAR0,SIGMA_TAU_0,SIGMA_NU_0,A,D0,D1,SIGMA_THETA,SIGMA_D,B,D_LOW,D_UP.
For the type "parametric", they are : 
N_ID_LOW,N_ID_HIGH,TMAX,DELTA,BY_LF,BY_HF,TAU_0,NU_0,SIGMA_OBS_HIGH,SIGMA_OBS_LOW,DMIN,DMAX,PAR0,A,D0,D1,SIGMA_THETA,SIGMA_D,B,D_LOW,D_UP.
Scripts for the simulation study with the current set_up are generated in folder "/domain/Rscripts".

The jobs can then be submitted with the bashscript "submit_jobs.sh". Before submitting the jobs, the current set_up may be checked by executing "check_current_settings.sh". Maximum 50 jobs can be submitted at a time. Project name, number of cores and nodes as well as walltime are set by changing the "submit_jobs.sh" script.
Once the calculation is done, the results of the parameter estimates are stored in .csv files in a folder "results_type_hyperparams" where hyperparams is the name of the hyperparameter file that was used for the simulation study, and type refers to the type of simulation study.
