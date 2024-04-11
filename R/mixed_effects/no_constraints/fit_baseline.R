

#SEED FOR REPRODUCTIBILITY
set.seed(42)

library(smoothSDE)
library(ggplot2)


########################### OUTPUT FILE TO KEEP RECORD OF THE SCRIPT EXECUTION #######################
filename=paste("fit_baseline_output.txt",sep="")

#Delete file if it exists
if (file.exists(filename)) {
  file.remove(filename)
}

file.create(filename)

####################################    GET NARWHAL DATA   ###########################################


# Set the path to the directory containing the data
par_dir=dirname(dirname(dirname(getwd()))) #parent directory 
narwhal_data_path <- file.path(par_dir,"Data", "Narwhals")

# DATA BEFORE EXPOSURE
dataBE1=read.csv(file.path(narwhal_data_path,"DataBE1.csv"), header = TRUE,dec = ".")

cat("Extracting trajectories before exposure and 1 day after tagging... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE1[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)


#  DATA BEFORE EXPOSURE WITH 12H  TAGGING EFFECT

dataBE2=read.csv(file.path(narwhal_data_path,"DataBE2.csv"), header = TRUE,dec = ".")

cat("Extracting all trajectories before exposure and 12 hours after tagging ... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE2[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)




##------------------------------------------------------------------------------------------------------



#Baseline models fitted on the data before the narwhals get in line of sight with the ship
# and 24 hours after they got tagged



##----------------------------------------------------------------------------------------------------

sigma_obs=0.045
n_obs=length(dataBE1$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))

######################## BASELINE MODEL WITHOUT RANDOM EFFECTS NOR COVARIATES #####################

#define model
formulas <- list(mu1 = ~1 ,mu2 =~1 ,tau = ~1,nu=~1)
par0 <- c(0,0,1,1)
baseline0 <- SDE$new(formulas = formulas,data = dataBE1,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))

#fit_model
baseline0$fit()
estimates0=as.list(baseline0$tmb_rep(),what="Est")
std0=as.list(baseline0$tmb_rep(),what="Std")

#############################   BASELINE MODEL WITH RANDOM EFFECTS IN ALL PAR   #############################

#define model
formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re") ,tau = ~s(ID,bs="re") ,nu=~s(ID,bs="re"))
par0 <- c(0,0,1,1)
baseline1 <- SDE$new(formulas = formulas,data = dataBE1,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))


#fit_model
baseline1$fit()
estimates1=as.list(baseline1$tmb_rep(),what="Est")
std1=as.list(baseline1$tmb_rep(),what="Std")

#plot parameters
res=baseline1$get_all_plots(baseline=NULL,model_name="baseline1",npost=1000,level=0.95)




#############################   BASELINE MODEL WITH RANDOM EFFECTS ONLY IN MU AND NU    #############################

#define model
formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re") ,tau = ~1 ,nu=~s(ID,bs="re"))
par0 <- c(0,0,1,1)
baseline2 <- SDE$new(formulas = formulas,data = dataBE1,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))


#fit_model
baseline2$fit()
estimates2=as.list(baseline2$tmb_rep(),what="Est")
std2=as.list(baseline2$tmb_rep(),what="Std")

#plot parameters
res=baseline2$get_all_plots(baseline=NULL,model_name="baseline2",npost=1000,level=0.95)


#############################   BASELINE MODEL WITH RANDOM EFFECTS ONLY IN NU AND FIXED MU    #############################

#define model
formulas <- list(mu1 = ~1,mu2 =~1 ,tau = ~1,nu=~s(ID,bs="re"))
par0 <- c(0,0,5,1)
baseline3 <- SDE$new(formulas = formulas,data = dataBE1,type = "CTCRW",response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),
                     other_data=list("log_sigma_obs0"=log(sigma_obs)))


#fit_model
baseline3$fit()
estimates3=as.list(baseline3$tmb_rep(),what="Est")
std3=as.list(baseline3$tmb_rep(),what="Std")

#plot parameters
res=baseline3$get_all_plots(baseline=NULL,model_name="baseline3",npost=1000,level=0.95)



#########################  RACVM BASELINE MODEL WITHOUT COVARIATES   ##############################

#define model
formulas <- list(mu1=~s(ID,bs="re"),mu2=~s(ID,bs="re"),
                 tau =~1,
                 nu=~s(ID,bs="re"),
                 omega=~s(ID,bs="re"))
par0 <- c(0,0,1,4,0)
baseline4<- SDE$new(formulas = formulas,data = dataBE1,type = "RACVM1",
                    response = c("x","y"),par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))



#fit_model
baseline4$fit()
estimates4=as.list(baseline4$tmb_rep(),what="Est")
std4=as.list(baseline4$tmb_rep(),what="Std")

res=baseline4$get_all_plots(baseline=NULL,model_name="baseline4",npost=1000,level=0.95)



###########################         AICs VALUES FOR THE BASELINE MODELS        ####################################

AIC_bas0=baseline0$AIC_marginal()
AIC_bas1=baseline1$AIC_marginal()
AIC_bas2=baseline2$AIC_marginal()
AIC_bas3=baseline3$AIC_marginal()
AIC_bas4=baseline4$AIC_marginal()

cat("BASELINE AICs VALUES WITH LOCATIONS OF FIRST DAY REMOVED TO AVOID TAGGING EFFECT",paste("without random effects :",AIC_bas0),
    paste("with random effects on tau, nu, mu :",AIC_bas1),paste("with random effects on nu, mu :",AIC_bas2),
    paste("with random effect on nu and mu=0:",AIC_bas3),
    paste("RACVM with random effects on mu and nu",AIC_bas4,"\n"),file=filename,sep="\n",append=TRUE)






##------------------------------------------------------------------------------------------------------



#Baseline models fitted on the data before the narwhals get in line of sight with the ship
# and 12 hours after they got tagged



##----------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------

sigma_obs=0.045
n_obs=length(dataBE2$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))

######################## BASELINE MODEL WITHOUT RANDOM EFFECTS NOR COVARIATES #####################

#define model
formulas <- list(mu1 = ~1 ,mu2 =~1 ,tau = ~1,nu=~1)
par0 <- c(0,0,1,4)
baseline0_12h <- SDE$new(formulas = formulas,data = dataBE2,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))

#fit_model
baseline0_12h$fit()
estimates0_12h=as.list(baseline0_12h$tmb_rep(),what="Est")
std0_12h=as.list(baseline0_12h$tmb_rep(),what="Std")

#############################   BASELINE MODEL WITH RANDOM EFFECTS IN ALL PAR   #############################

#define model
formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re") ,tau = ~s(ID,bs="re") ,nu=~s(ID,bs="re"))
par0 <- c(0,0,1,4)
baseline1_12h <- SDE$new(formulas = formulas,data = dataBE2,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))


#fit_model
baseline1_12h$fit()
estimates1_12h=as.list(baseline1_12h$tmb_rep(),what="Est")
std1_12h=as.list(baseline1_12h$tmb_rep(),what="Std")

#plot parameters
res=baseline1_12h$get_all_plots(baseline=NULL,model_name="baseline1_12h",npost=1000,level=0.95)




#############################   BASELINE MODEL WITH RANDOM EFFECTS ONLY IN MU AND NU    #############################

#define model
formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re") ,tau = ~1 ,nu=~s(ID,bs="re"))
par0 <- c(0,0,1,4)
baseline2_12h <- SDE$new(formulas = formulas,data = dataBE2,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))


#fit_model
baseline2_12h$fit()
estimates2_12h=as.list(baseline2_12h$tmb_rep(),what="Est")
std2_12h=as.list(baseline2_12h$tmb_rep(),what="Std")

#plot parameters
res=baseline2_12h$get_all_plots(baseline=NULL,model_name="baseline2_12h",npost=1000,level=0.95)


#############################   BASELINE MODEL WITH RANDOM EFFECTS ONLY IN NU AND FIXED MU    #############################

#define model
formulas <- list(mu1 = ~1,mu2 =~1 ,tau = ~1,nu=~s(ID,bs="re"))
par0 <- c(0,0,1,4)
baseline3_12h<- SDE$new(formulas = formulas,data = dataBE2,type = "CTCRW",response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),
                        other_data=list("log_sigma_obs0"=log(sigma_obs)))


#fit_model
baseline3_12h$fit()
estimates3_12h=as.list(baseline3_12h$tmb_rep(),what="Est")
std3_12h=as.list(baseline3_12h$tmb_rep(),what="Std")

#plot parameters
res=baseline3_12h$get_all_plots(baseline=NULL,model_name="baseline3_12h",npost=1000,level=0.95)



#########################  RACVM BASELINE MODEL WITHOUT COVARIATES   ##############################

#define model
formulas <- list(mu1=~s(ID,bs="re"),mu2=~s(ID,bs="re"),
                 tau =~1,
                 nu=~s(ID,bs="re"),
                 omega=~s(ID,bs="re"))
par0 <- c(0,0,1,4,0)
baseline4_12h<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM1",
                    response = c("x","y"),par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))



#fit_model
baseline4_12h$fit()
estimates4_12h=as.list(baseline4_12h$tmb_rep(),what="Est")
std4_12h=as.list(baseline4_12h$tmb_rep(),what="Std")

res=baseline4_12h$get_all_plots(baseline=NULL,model_name="baseline4_12h",npost=1000,level=0.95)


###########################         AICs VALUES FOR THE BASELINE MODELS        ####################################

AIC_bas0_12h=baseline0_12h$AIC_marginal()
AIC_bas1_12h=baseline1_12h$AIC_marginal()
AIC_bas2_12h=baseline2_12h$AIC_marginal()
AIC_bas3_12h=baseline3_12h$AIC_marginal()
AIC_bas4_12h=baseline4_12h$AIC_marginal()

cat("BASELINE AICs VALUES WITH LOCATIONS OF FIRST 12h REMOVED TO AVOID TAGGING EFFECT",paste("without random effects :",AIC_bas0_12h),
    paste("with random effects on tau, nu, mu :",AIC_bas1_12h),paste("with random effects on nu, mu :",AIC_bas2_12h),
    paste("with random effect on nu and mu=0:",AIC_bas3_12h),
    paste("RACVM with random effects on mu and nu",AIC_bas4_12h),file=filename,sep="\n",append=TRUE)



