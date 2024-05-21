# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-17
#
# Script Name:    ~/Documents/Research/Projects/narwhals_smoothSDE/R/application/CRCVM/fit_response.R
#
# Script Description: Fit response CRCVM models
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console

#SEED FOR REPRODUCTIBILITY
set.seed(42)

library(smoothSDE)
library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)




################################### OUTPUT FILE TO KEEP RECORD OF THE SCRIPT EXECUTION ####################
filename=paste("response_output.txt",sep="")

#Delete file if it exists
if (file.exists(filename)) {
  file.remove(filename)
}

file.create(filename)

########################################    GET NARWHAL DATA   #############################################


# Set the path to the directory containing the data
par_dir=dirname(dirname(dirname(getwd()))) #parent directory 
narwhal_data_path <- file.path(par_dir,"Data", "Narwhals")  


# DATA AFTER EXPOSURE

dataAE=read.csv(file.path(narwhal_data_path,"DataAE.csv"), header = TRUE,dec = ".")


cat("Extracting trajectories after exposure.. ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataAE[,1]),"positions measured after exposure \n"),file=filename,sep="\n",append=TRUE)


############################## SOME PREPROCESSING FOR EXPSHORE ####################################

D_low=0.07
D_up=3
dataAE[dataAE$DistanceShore> D_up,"ExpShore"]=0
dataAE[dataAE$DistanceShore< D_low,"ExpShore"]=1/D_low

n_obs=length(dataAE$time)



########### MODEL WITH ONLY ANGLE NORMAL ##################
par0 <- c(0,0,1,4,0)

sigma_obs=0.05
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))

#model formula
formulas <- list(mu1=~1,mu2=~1,
                 tau =~s(ID,bs="re"),
                 nu=~s(ID,bs="re"),
                 omega=~s(AngleNormal,k=4,bs="cs"))

response1<- SDE$new(formulas = formulas,data = dataAE,type = "RACVM",response = c("x","y"),
                    par0 = par0,other_data=list("log_sigma_obs0"=log(0.05)),fixpar=c("mu1","mu2"))

#fit_model
response1$fit()
estimates1=as.list(response1$tmb_rep(),what="Est")
std1=as.list(response1$tmb_rep(),what="Std")

res=response1$get_all_plots(baseline=NULL,model_name="response1",show_CI="pointwise",save=TRUE)

#we see a clear effect of the angle on both omega and tau

############## MODEL WITH FIXED SPLINE OF ANGLE NORMAL AND TAU AND NU INTERCEPTS ##################

par0 <- c(0,0,1,4,0)

sigma_obs=0.05
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))

#model formula
formulas <- list(mu1=~1,mu2=~1,
                 tau =~s(ExpShip,k=3,bs="tp")+s(ID,bs="re"),
                 nu=~s(ExpShip,k=3,bs="tp")+s(ID,bs="re"),
                 omega=~s(AngleNormal,k=4,bs="cs"))

response2<- SDE$new(formulas = formulas,data = dataAE,type = "RACVM",
                    response = c("x","y"),par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),
                    fixpar=c("mu1","mu2"),map=list(coeff_re=factor(c(1:2,rep(NA,6),3:4,rep(NA,6),rep(NA,3))),coeff_fe=factor(rep(NA,5))))

new_coeff_re=c(rep(0,2),baseline1$coeff_re()[paste("tau.s(ID).",1:6,sep=""),1],rep(0,2),baseline1$coeff_re()[paste("nu.s(ID).",1:6,sep=""),1],
              baseline1$coeff_re()[paste("omega.s(AngleNormal).",1:3,sep=""),1])

response2$update_coeff_re(new_coeff_re)
response2$update_coeff_fe(baseline1$coeff_fe()[,1])


#fit_model
response2$fit()
estimates_res2=as.list(response2$tmb_rep(),what="Est")
std_res2=as.list(response2$tmb_rep(),what="Std")


#plot parameters

#plot parameters
xmin=list("AngleNormal"=-pi,"ExpShip"=1/45)
xmax=list("AngleNormal"=pi,"ExpShip"=1/3)
link=list("ExpShip"=(\(x) 1/x))
xlabel=list("ExpShip"="Distance to ship")

res=response2$get_all_plots(baseline=baseline1,model_name="response2",xmin=xmin,
                            xmax=xmax,link=link,xlabel=xlabel,show_CI="none",save=TRUE)

####### FULL MODEL WITH FIXED BASELINE SPLINE COEFFS

#define model
formulas <- list(mu1=~1,mu2=~1,tau =~s(ExpShip,k=3,bs="cs")+s(ID,bs="re"),
                 nu=~s(ExpShip,k=3,bs="cs")+s(ID,bs="re"),
                 omega=~ti(ExpShore,k=5,bs="cs")+ti(AngleNormal,k=5,bs="cs")+ti(ExpShore,AngleNormal,k=5,bs="cs")+s(ID,bs="re"))
par0 <- c(0,0,1,4,0)
response4<- SDE$new(formulas = formulas,data = dataAE,type = "RACVM",
                    response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),other_data=list("log_sigma_obs0"=log(sigma_obs)),
                    map=list(coeff_re=factor(c(1:2,rep(NA,6),3:4,rep(NA,6),rep(NA,24),rep(NA,6))),coeff_fe=factor(rep(NA,5))))


new_coeff_re=c(rep(0,2),baseline2$coeff_re()[paste("tau.s(ID).",1:6,sep=""),1],rep(0,2),baseline2$coeff_re()[paste("nu.s(ID).",1:6,sep=""),1],
               baseline2$coeff_re()[paste("omega.ti(ExpShore).",1:4,sep=""),1],baseline2$coeff_re()[paste("omega.ti(AngleNormal).",1:4,sep=""),1],
          baseline2$coeff_re()[paste("omega.ti(ExpShore,AngleNormal).",1:16,sep=""),1],baseline2$coeff_re()[paste("omega.s(ID).",1:6,sep=""),1])

response4$update_coeff_re(new_coeff_re)
response4$update_coeff_fe(baseline2$coeff_fe()[,1])


#fit_model
response4$fit()
estimates4=as.list(response4$tmb_rep(),what="Est")
std4=as.list(response4$tmb_rep(),what="Std")

#plot parameters
xmin=list("ExpShore"=1/D_up,"AngleNormal"=-pi,"ExpShip"=1/45)
xmax=list("ExpShore"=1/D_low,"AngleNormal"=pi,"ExpShip"=1/3)
link=list("ExpShore"=(\(x) 1/x),"ExpShip"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore","ExpShip"="Distance to ship")

response4$get_all_plots(baseline=baseline2,model_name="response4",xmin=xmin,
                        xmax=xmax,link=link,xlabel=xlabel,show_CI="none",save=TRUE)


#############################      RESPONSE AIC VALUES    #####################################################

AIC_resp1=response1$AIC_marginal()
AIC_resp2=response2$AIC_marginal()
AIC_resp3=response3$AIC_marginal()




#################################    POSTERIOR PREDICTIVE CHECKS  #########################################


check_response2=response2$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response2,"response2")

check_response3=response3$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response3,"response3")

check_response4=response4$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response4,"response4")

check_response5=response5$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response5,"response5")

check_response6=response6$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response6,"response6")

check_response7=response7$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response7,"response7")

check_response8=response8$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response8,"response8")

check_response9=response9$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response9,"response9")

check_response10=response10$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response10,"response10")

check_response11=response11$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response11,"response11")

check_response12=response12$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response12,"response12")

check_response13=response13$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response13,"response13")





check_response2offset=response2offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response2offset,"response2offset")

check_response3offset=response3offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response3offset,"response3offset")

check_response4offset=response4offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response4offset,"response4offset")

check_response5offset=response5offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response5offset,"response5offset")

check_response6offset=response6offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response6,"response6offset")

check_response7offset=response7offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response7offset,"response7offset")

check_response8offset=response8offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response8offset,"response8offset")

check_response9offset=response9offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response9offset,"response9offset")

check_response10offset=response10offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response10offset,"response10offset")

check_response11offset=response11offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response11offset,"response11offset")

check_response12offset=response12offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response12offset,"response12offset")

check_response13offset=response13offste$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response13offset,"response13offset")







#######################   SIMULATE DATA WITH FITTED PARAMETERS MODEL AND  TRY TO RECOVER IT #########################

# models without offset
check_fe_estimations(ctcrw=response2,model_name="response2",links=list("ExpShip"=(\(x) 1/x)),
                      xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response3,model_name="response3",links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response4,model_name="response4",links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)





#models with offset
check_fe_estimations(ctcrw=response2offset,model_name="response2offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response3offset,model_name="response3offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response4offset,model_name="response4offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)


check_fe_estimations(ctcrw=response8offset,model_name="response8offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response9offset,model_name="response9offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)




