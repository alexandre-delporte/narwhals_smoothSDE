



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


# DATA BEFORE EXPOSURE

dataBE1=read.csv(file.path(narwhal_data_path,"DataBE1.csv"), header = TRUE,dec = ".")

cat("Extracting trajectories before exposure and 1 day after tagging... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE1[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)

dataBE2=read.csv(file.path(narwhal_data_path,"DataBE2.csv"), header = TRUE,dec = ".")



cat("Extracting trajectories before exposure and 12h after tagging... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE2[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)


# DATA AFTER EXPOSURE

dataAE=read.csv(file.path(narwhal_data_path,"DataAE.csv"), header = TRUE,dec = ".")


cat("Extracting trajectories after exposure.. ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataAE[,1]),"positions measured after exposure \n"),file=filename,sep="\n",append=TRUE)



############################## SOME PREPROCESSING FOR EXPSHORE ####################################

#keep points far enough from land
dataAE=dataAE[dataAE$DistanceShore>0.05,]

#put threshold on exposure
dataAE[dataAE$DistanceShore>1,"ExpShore"]=0


############################################# SET MEASUREMENT ERROR ###############################################

sigma_obs=0.045
n_obs=length(dataAE$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))


#######################  RACVM  MODEL WITH tensor splines te of AngleNormal and Distance Shore in omega #######################

#define model
formulas <- list(mu1=~1,mu2=~1,
                 tau =~1,
                 nu=~1,omega=~te(ExpShore,AngleNormal,k=5,bs="cs"))
par0 <- c(0,0,1,4,0)
response1<- SDE$new(formulas = formulas,data = dataAE,type = "RACVM1",
                    response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),
                    other_data=list("H"=H))



#fit_model
response1$fit()
estimates1=as.list(response1$tmb_rep(),what="Est")
std1=as.list(response1$tmb_rep(),what="Std")


#plot parameters
xmin=list("ExpShore"=1/1,"AngleNormal"=-pi+pi/20)
xmax=list("ExpShore"=1/0.05,"AngleNormal"=pi-pi/20)
link=list("ExpShore"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore")


res=response1$get_all_plots(baseline=NULL,model_name="response1",xmin=xmin,
              xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)

p_far=response1$plot_par(var="AngleNormal",par_names=c("omega"),covs=data.frame("ExpShore"=1/0.5,"AngleNormal"=0),show_CI="pointwise")
p_close=response1$plot_par(var="AngleNormal",par_names=c("omega"),covs=data.frame("ExpShore"=1/0.1,"AngleNormal"=0),show_CI="pointwise")


#########################  RACVM  MODEL WITH tensor splines of  AngleNormal, Distance shore in omega ##############################

#define model
formulas <- list(mu1=~1,mu2=~1,tau =~s(ExpShore,k=5,bs="cs")+s(AngleNormal,k=5,bs="cs"),
                 nu=~s(ID,bs="re"),
                 omega=~ti(AngleNormal,k=5,bs="cs")+ti(ExpShore,k=5,bs="cs")+ti(AngleNormal,ExpShore,k=5,bs="cs"))
par0 <- c(0,0,1,1,0)
response2<- SDE$new(formulas = formulas,data =s dataAE,type = "RACVM1",
                    response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),other_data=list("H"=H))



#fit_model
response2$fit()
estimates2=as.list(response2$tmb_rep(),what="Est")
std2=as.list(response2$tmb_rep(),what="Std")


#plot parameters

response2$get_all_plots(baseline=NULL,model_name="response2",xmin=xmin,
                        xmax=xmax,xlabel=xlabel,npost=1000,level=0.95,link=link)

p_far=response2$plot_par(var="AngleNormal",par_names=c("omega"),covs=data.frame("ExpShore"=1/0.5,"AngleNormal"=0),show_CI="pointwise")
p_close=response2$plot_par(var="AngleNormal",par_names=c("omega"),covs=data.frame("ExpShore"=1/0.1,"AngleNormal"=0),show_CI="pointwise")
#########################  RACVM  MODEL WITH additive splines of  AngleNormal and ExpShore in omega and tau ##############################

#define model
formulas <- list(mu1=~1,mu2=~1,tau =~1,
                 nu=~s(ExpShip,k=10,bs="cs"),
                 omega=~ti(AngleNormal,k=5,bs="cs")+ti(ExpShore,k=5,bs="cs")+ti(AngleNormal,ExpShore,k=5,bs="cs"))
par0 <- c(0,0,1,4,0)
response3<- SDE$new(formulas = formulas,data = dataAE,type = "RACVM1",
                    response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),other_data=list("H"=H))



#fit_model
response3$fit()
estimates3=as.list(response3$tmb_rep(),what="Est")
std3=as.list(response3$tmb_rep(),what="Std")

#plot parameters
xmin=list("ExpShore"=1/0.5,"AngleNormal"=-pi+pi/20,"ExpShip"=1/80)
xmax=list("ExpShore"=1/0.05,"AngleNormal"=pi-pi/20,"ExpShip"=1/0.5)
link=list("ExpShore"=(\(x) 1/x),"ExpShip"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore","ExpShip"="Distance to ship")

response3$get_all_plots(baseline=NULL,model_name="response3",xmin=xmin,
              xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)



#########################  RACVM  MODEL WITH additive splines of  AngleNormal and ExpShore in omega and tau ##############################

#define model
formulas <- list(mu1=~1,mu2=~1,tau =~s(ExpŜhip,k=10,bs="cs"),
                 nu=~s(ExpŜhip,k=10,bs="cs")+s(ID,bs="re"),
                 omega=~ti(AngleNormal,k=5,bs="cs")+ti(ExpShore,k=5,bs="cs")+ti(AngleNormal,ExpShore,k=5,bs="cs"))
par0 <- c(0,0,1,1,0)
response4<- SDE$new(formulas = formulas,data = dataAE,type = "RACVM1",
                    response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),other_data=list("H"=H))



#fit_model
response4$fit()
estimates4=as.list(response4$tmb_rep(),what="Est")
std4=as.list(response4$tmb_rep(),what="Std")

#plot parameters


response4$get_all_plots(baseline=NULL,model_name="response4",xmin=xmin,
                        xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)


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




