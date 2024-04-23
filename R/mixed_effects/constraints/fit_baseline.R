

#SEED FOR REPRODUCTIBILITY
set.seed(42)

library(smoothSDE)
library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)





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




#Baseline models fitted on the data before the narwhals get in line of sight with the ship
# and 24 hours after they got tagged

############################## SOME PREPROCESSING FOR EXPSHORE ####################################

#keep points far enough from land
dataBE2=dataBE2[dataBE2$DistanceShore>0.03,]
dataBE1=dataBE1[dataBE1$DistanceShore>0.03,]

#put threshold on exposure
dataBE1[dataBE1$DistanceShore>5,"ExpShore"]=0
dataBE2[dataBE2$DistanceShore>5,"ExpShore"]=0


############################################# SET MEASUREMENT ERROR ###############################################

sigma_obs=0.03
n_obs=length(dataBE2$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))

#########################  ANGLE NORMAL IN OMEGA   ##############################

#initial parameters
par0 <- c(0,0,1,4,0)

#model formula
formulas <- list(mu1=~1,mu2=~1,
                 tau =~1,
                 nu=~s(ID,bs="re"),
                 omega=~s(AngleNormal,k=10,bs="cs"))

baseline1<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM1",response = c("x","y"),
                    par0 = par0,other_data=list("H"=H),fixpar=c("mu1","mu2"))

#fit_model
baseline1$fit()
estimates1=as.list(baseline1$tmb_rep(),what="Est")
std1=as.list(baseline1$tmb_rep(),what="Std")

res=baseline1$get_all_plots(baseline=NULL,model_name="baseline1",npost=1000,level=0.95,show_CI=TRUE)


#########################  TE SPLINES ANGLE NORMAL AND EXP SHORE IN OMEGA   ##############################

formulas <- list(mu1 = ~1 ,mu2 =~1,tau = ~1,
                 nu=~s(ID,bs="re"),omega=~te(AngleNormal,ExpShore,k=c(6,3),bs="cs"))

baseline2<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM1",
                    response = c("x","y"),par0 = par0,other_data=list("H"=H),
                    fixpar=c("mu1","mu2"))

#fit_model
baseline2$fit()
estimates2=as.list(baseline2$tmb_rep(),what="Est")
std2=as.list(baseline2$tmb_rep(),what="Std")


#plot parameters
xmin=list("ExpShore"=1/5,"AngleNormal"=-pi+pi/20)
xmax=list("ExpShore"=1/0.05,"AngleNormal"=pi-pi/20)
link=list("ExpShore"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore")

res=baseline2$get_all_plots(baseline=NULL,model_name="baseline2",npost=1000,level=0.95,
                        xmin=xmin,xmax=xmax,link=link,xlabel=xlabel)

p_close=baseline2$plot_par(var="AngleNormal",par_names=c("omega"),covs=data.frame("ExpShore"=0.5,"AngleNormal"=0),show_CI="pointwise")

#########################  BIVARIATE SPLINES ANGLE NORMAL AND EXPSHORE IN OMEGA  ##############################

formulas <- list(mu1 = ~1 ,mu2 =~1,tau = ~1,
                 nu=~1,omega=~ti(AngleNormal,k=5,bs="cs")+ti(ExpShore,k=5,bs="cs")+
                   ti(AngleNormal,ExpShore,k=5,bs="cs"))

baseline3<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM1",
                    response = c("x","y"),par0 = par0,other_data=list("log_sigma_obs0"=log(0.04)),
                    fixpar=c("mu1","mu2"))

#fit_model
baseline3$fit()
estimates3=as.list(baseline3$tmb_rep(),what="Est")
std3=as.list(baseline3$tmb_rep(),what="Std")


res=baseline3$get_all_plots(baseline=NULL,model_name="baseline3",npost=1000,level=0.95,
                        xmin=xmin,xmax=xmax,link=link,xlabel=xlabel)

p_close=baseline3$plot_par(var="AngleNormal",par_names=c("omega"),covs=data.frame("ExpShore"=0.5,"AngleNormal"=0),show_CI="pointwise")

###########################         AICs VALUES FOR THE BASELINE MODELS        ####################################

AIC_bas1=baseline1$AIC_marginal()
AIC_bas2=baseline2$AIC_marginal()
AIC_bas3=baseline3$AIC_marginal()



cat("--------------------- BASELINE MODELS ---------------------: \n",
    baseline1$formulas(),AIC_bas1,
    baseline2$formulas(),AIC_bas2,baseline3$formulas(),AIC_bas3,file=filename,sep="\n",append=TRUE)
