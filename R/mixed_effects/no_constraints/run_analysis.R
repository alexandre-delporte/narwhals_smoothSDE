
cat("\n \n ***** FIT BASELINE MODELS ***** \n \n")

source("fit_baseline.R")

cat("\n \n ***** CHECK BASELINE MODELS  ***** \n \n")

source("check_baseline_models.R")

cat("\n \n ***** SIMULATE FROM FITTED BASELINE MODELS ***** \n \n")
source("simulate_baseline_models.R")

cat("\n \n ***** FIT RESPONSE MODELS ***** \n \n")
source("fit_response_baseline2.R")

source("fit_response_baseline3.R")

cat("\n \n ***** CHECK RESPONSE MODELS ***** \n \n")
source("check_response_models.R")

cat("***** \n \n SIMULATE FROM FITTED RESPONSE MODELS ***** \n \n")
source("simulate_response_models.R")