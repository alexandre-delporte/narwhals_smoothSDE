

#SEED FOR REPRODUCTIBILITY
set.seed(42)

#REQUIRED LIBRARIES
library(smoothSDE)
library(ggplot2)
#source(file.path("/home", "delporta","Documents","Recherche","Codes","smoothSDE","R","utility.R"))

########################################    MANAGE DIRECTORIES  #################################

#directory of mixed effects analysis
path <- file.path("/home", "delporta","Documents","Recherche","Codes",
                  "narwhals_movement_analysis","smoothSDE","mixed_effects","no_constraints")
setwd(path)


################################### OUTPUT FILE TO KEEP RECORD OF THE SCRIPT EXECUTION ####################
filename=paste("fit_response_output.txt",sep="")

#Delete file if it exists
if (file.exists(filename)) {
  file.remove(filename)
}

file.create(filename)

########################################    GET NARWHAL DATA   #############################################

data_path<- file.path("/home","delporta","Documents","Recherche","Données","narvals")       

# DATA BEFORE EXPOSURE

dataBE1=read.csv(file.path(data_path,"DataBE1.csv"), header = TRUE,dec = ".")

cat("Extracting trajectories before exposure and 1 day after tagging... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE1[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)

dataBE2=read.csv(file.path(data_path,"DataBE2.csv"), header = TRUE,dec = ".")

cat("Extracting trajectories before exposure and 12h after tagging... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE2[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)


# DATA AFTER EXPOSURE

dataAE=read.csv(file.path(data_path,"DataAE.csv"), header = TRUE,dec = ".")

cat("Extracting trajectories after exposure.. ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataAE[,1]),"positions measured after exposure \n"),file=filename,sep="\n",append=TRUE)




############################################# SET MEASUREMENT ERROR ###############################################

sigma_obs=0.03
n_obs=length(dataAE$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))

####################################  RESPONSE MODELS WITH ExpShip IN NU ###########################################


#initialize parameters with estimated baseline3 values
par0=sapply(1:4, function(i) {baseline3$invlink()[[i]](estimates3$coeff_fe[i])})

#Find best degree of freedom for cubic splines

list_k=c(3,4,7,10,15,20,30)
AIC_df7=data.frame("k"=list_k,"AIC_marginal"=rep(NA,length(list_k)),"AIC_conditional"=rep(NA,length(list_k)))

for (i in 1:length(list_k)) {
  
  formulas <- list(mu1 = ~1 ,mu2 =~1,tau = ~1,
                   nu=~s(ExpShip,k=list_k[i],bs="cs")+s(ID,bs="re"))

  try_response7 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                       par0 = par0,other_data=list("H"=H),fixpar=c("mu1","mu2"))
  try_response7$fit()
  
  AIC_df7[i,"AIC_marginal"]=try_response7$AIC_marginal()
  AIC_df7[i,"AIC_conditional"]=try_response7$AIC_conditional()
}

#plot AIC as a function of number of spline knots
plotAIC7=ggplot(data=AIC_df7,aes(x=k))+geom_line(aes(y=AIC_marginal,color="AIC_marginal"))+geom_point(aes(y=AIC_marginal),shape=21)+
  geom_line(aes(k, AIC_conditional,color="AIC_conditional"))+geom_point(aes(k, AIC_conditional),shape=21)+
  xlab("Cubic splines degree of freedom")+ylab("AIC")+
  scale_color_manual(values = c("AIC_marginal" = "red", "AIC_conditional" = "blue"),
                     labels = c("Marginal AIC", "Conditional AIC"),name=NULL)
  
ggsave(path="response7","plotAIC7.png",plot=plotAIC7)




#set k to best degree of freedom
formulas <- list(mu1 = ~1 ,mu2 =~1,tau = ~1,
                 nu=~s(ExpShip,k=6,bs="cs")+s(ID,bs="re"))

#FIT MODEL WITHOUT OFFSET

response7 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),fixpar=c("mu1","mu2"))
response7$fit()


xmin=list("ExpShip"=1/80)
xmax=list("ExpShip"=1/3)
link=list("ExpShip"=(\(x) 1/x))
xlabel=list("ExpShip"="Distance to ship")

plot_resp7=response1$get_all_plots(baseline=baseline3,model_name="response7",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95,show_CI=TRUE)


## MODEL WITH OFFSET AS IN BASELINE 3

map=list("coeff_fe"=factor(rep(NA,4)))

response7offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",
                           response = c("x","y"),par0 = par0,map=map,other_data=list("H"=H),
                           fixpar=c("mu1","mu2"))
response7offset$fit()

plot_resp7offset=response7offset$get_all_plots(baseline=baseline3,model_name="response7offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)


##################################  RESPONSE MODELS WITH ExpShip IN TAU  #################################################


#Find best degree of freedom for cubic splines

AIC_df8=data.frame("k"=list_k,"AIC_marginal"=rep(NA,length(list_k)),"AIC_conditional"=rep(NA,length(list_k)))

for (i in 1:length(list_k)) {

  formulas <- list(mu1 = ~1 ,mu2 =~1 ,tau =~s(ExpShip,k=list_k[i],bs="cs"),
                   nu=~s(ID,bs="re"))
  
  try_response8 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,other_data=list("H"=H),fixpar=c("mu1","mu2"))
  try_response8$fit()
  
  AIC_df8[i,"AIC_marginal"]=try_response8$AIC_marginal()
  AIC_df8[i,"AIC_conditional"]=try_response8$AIC_conditional()
}

#PLot AIC as a function of number of spline knots
plotAIC8=ggplot(data=AIC_df8,aes(x=k))+geom_line(aes(y=AIC_marginal,color="AIC_marginal"))+geom_point(aes(y=AIC_marginal),shape=21)+
  geom_line(aes(k, AIC_conditional,color="AIC_conditional"))+geom_point(aes(k, AIC_conditional),shape=21)+
  xlab("Cubic splines degree of freedom")+ylab("AIC")+
  scale_color_manual(values = c("AIC_marginal" = "red", "AIC_conditional" = "blue"),
                     labels = c("Marginal AIC", "Conditional AIC"),name=NULL)

ggsave(path="response8","plotAIC8.png",plot=plotAIC8)


#set k to best degree of freedom
formulas <- list(mu1 = ~1 ,mu2 =~1,tau =~s(ExpShip,k=7,bs="cs"),
                 nu=~s(ID,bs="re"))

## FIT MODEL WITHOUT OFFSET

response8 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),fixpar=c("mu1","mu2"))
response8$fit()

#get parameters plots
plot_resp8=response8$get_all_plots(baseline=baseline3,model_name="response8",xmin=xmin,xmax=xmax,
              link=link,xlabel=xlabel,npost=1000,level=0.95)

## FIT MODEL WITH OFFSET AS IN BASELINE 3
response8offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,other_data=list("H"=H),map=map,fixpar=c("mu1","mu2"))
response8offset$fit()

#get parameters plot
plot_resp8offset=response8offset$get_all_plots(baseline=baseline3,model_name="response8offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)


###################################    RESPONSE MODELS WITH ExpShip IN TAU AND NU ##############################################



#Find best degree of freedom for cubic splines

AIC_df9=data.frame("k"=list_k,"AIC_marginal"=rep(NA,length(list_k)),"AIC_conditional"=rep(NA,length(list_k)))

for (i in 1:length(list_k)) {
  
  
  formulas <- list(mu1 = ~1 ,mu2 =~1 ,tau = ~s(ExpShip,k=list_k[i],bs="cs"),
                   nu=~s(ExpShip,k=list_k[i],bs="cs")+s(ID,bs="re"))
  
  try_response9 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,other_data=list("H"=H),fixpar=c("mu1","mu2"))
  try_response9$fit()
  
  AIC_df9[i,"AIC_marginal"]=try_response9$AIC_marginal()
  AIC_df9[i,"AIC_conditional"]=try_response9$AIC_conditional()
}

plotAIC9=ggplot(data=AIC_df9,aes(x=2*k))+geom_line(aes(y=AIC_marginal,color="AIC_marginal"))+geom_point(aes(y=AIC_marginal),shape=21)+
  geom_line(aes(y=AIC_conditional,color="AIC_conditional"))+geom_point(aes(y=AIC_conditional),shape=21)+
  xlab("Number of coefficients to estimate (=2q)")+ylab("AIC")+
  scale_color_manual(values = c("AIC_marginal" = "red", "AIC_conditional" = "blue"),
                     labels = c("Marginal AIC", "Conditional AIC"),name=NULL)

ggsave(path="response9","plotAIC9.png",plot=plotAIC9,width=10,height=5)


#Set k to best degree of freedom
formulas <- list(mu1 = ~1 ,mu2 =~1,tau = ~s(ExpShip,k=10,bs="cs"),
                 nu=~s(ExpShip,k=10,bs="cs")+s(ID,bs="re"))



## FIT MODEL WITHOUT OFFSET
response9 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),fixpar=c("mu1","mu2"))
response9$fit()

plot_resp9=response9$get_all_plots(baseline=baseline3,model_name="response9",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95,show_CI=TRUE)

#inspect specific values of tau
tau_data=plot_resp9$data_tau_ExpShip
q1=quantile(dataAE$Dist_Ship,0.05,na.rm=TRUE)
q2=quantile(dataAE$Dist_Ship,0.25,na.rm=TRUE)
q3=quantile(dataAE$Dist_Ship,0.75,na.rm=TRUE)
very_close_tau=mean(tau_data[1/tau_data$ExpShip<q1,"par"])
close_tau=mean(tau_data[1/tau_data$ExpShip>q1 & 1/tau_data$ExpShip<q2,"par"])
medium_tau=mean(tau_data[1/tau_data$ExpShip>q2 & 1/tau_data$ExpShip<q3,"par"])
far_tau=mean(tau_data[1/tau_data$ExpShip>q3,"par"])

#inspect specific values of nu
nu_data=plot_resp9$data_nu_ExpShip
very_close_nu=mean(nu_data[1/nu_data$ExpShip<q1,"par"])
close_nu=mean(nu_data[1/nu_data$ExpShip>q1 & 1/nu_data$ExpShip<q2,"par"])
medium_nu=mean(nu_data[1/nu_data$ExpShip>q2 & 1/nu_data$ExpShip<q3,"par"])
far_nu=mean(nu_data[1/nu_data$ExpShip>q3,"par"])

## FIT MODEL WITH OFFSET AS IN BASELINE 2
response9offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,map=map,other_data=list("H"=H),fixpar=c("mu1","mu2"))
response9offset$fit()

plot_resp9offset=response9offset$get_all_plots(baseline=baseline3,model_name="response9offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)


###############################   RESPONSE MODELS WITH XT AND XI COVARIATE IN NU #############################################

formulas <- list(mu1 = ~1 ,mu2 =~1,tau = ~1,
                 nu=~s(ExpShipT,k=7,bs="cs")+
                   s(ExpShipI,k=7,bs="cs")+s(ID,bs="re"))

## MODEL WITHOUT OFFSET
response10 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),fixpar=c("mu1","mu2"))
response10$fit()


xmin=list("ExpShipT"=1/80,"ExpShipI"=1/80)
xmax=list("ExpShipT"=1/3,"ExpShipI"=1/3)
link=list("ExpShipT"=(\(x) 1/x),"ExpShipI"=(\(x) 1/x))
xlabel=list("ExpShipT"="Distance to ship","ExpShipI"="Distance to ship")

plot_resp10=response10$get_all_plots(baseline=baseline3,model_name="response10",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)


## MODEL WITH OFFSET AS IN BASELINE 1

response10offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,map=map,other_data=list("H"=H),fixpar=c("mu1","mu2"))
response10offset$fit()

plot_resp10offset=response10offset$get_all_plots(baseline=baseline3,model_name="response10offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)


###################################  RESPONSE MODEL WITH ExpShipT AND ExpShipI COVARIATES IN TAU  ####################################


formulas <- list(mu1 = ~1 ,mu2 =~1 ,tau =~s(ExpShipT,k=7,bs="cs")+s(ExpShipI,k=7,bs="cs"),
                 nu=~s(ID,bs="re"))

## MODEL WITHOUT OFFSET
response11 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),fixpar=c("mu1","mu2"))
response11$fit()


plot_resp11=response11$get_all_plots(baseline=baseline3,model_name="response5",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)


## MODEL WITH OFFSET AS IN BASELINE 1

response11offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,other_data=list("H"=H),map=map,fixpar=c("mu1","mu2"))
response11offset$fit()


plot_resp11offset=response11offset$get_all_plots(baseline=baseline3,model_name="response11offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)



#############################   RESPONSE  MODEL WITH ExpShipT AND ExpShipI COVARIATES IN NU AND TAU    ####################################


formulas <- list(mu1 = ~1 ,mu2 =~1 ,tau = ~s(ExpShipT,k=7,bs="cs")+s(ExpShipI,k=7,bs="cs")
                 ,nu=~s(ExpShipT,k=7,bs="cs")+s(ExpShipI,k=7,bs="cs")+s(ID,bs="re"))

## MODEL WITHOUT OFFSET

response12 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("H"=H),fixpar=c("mu1","mu2"))
response12$fit()


plot_resp12=response12$get_all_plots(baseline=baseline3,model_name="response12",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)


## MODEL WITH OFFSET AS IN BASELINE 1

response12offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,map=map,other_data=list("H"=H),fixpar=c("mu1","mu2"))
response12offset$fit()


plot_resp12offset=response12offset$get_all_plots(baseline=baseline3,model_name="response12offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,npost=1000,level=0.95)


#############################      RESPONSE AIC VALUES    #####################################################

AIC_resp7=response7$AIC_marginal()
AIC_resp8=response8$AIC_marginal()
AIC_resp9=response9$AIC_marginal()
AIC_resp10=response10$AIC_marginal()
AIC_resp11=response11$AIC_marginal()
AIC_resp12=response12$AIC_marginal()



AIC_resp7offset=response7offset$AIC_marginal()
AIC_resp8offset=response8offset$AIC_marginal()
AIC_resp9offset=response9offset$AIC_marginal()
AIC_resp10offset=response10offset$AIC_marginal()
AIC_resp11offset=response11offset$AIC_marginal()
AIC_resp12offset=response12offset$AIC_marginal()



cat("--------------------- RESPONSE MODELS BASED ON BASELINE 3 ---------------------: \n",
    "WITHOUT OFFSET : \n",response7$message(),AIC_resp7,
    response8$message(),AIC_resp8,response9$message(),AIC_resp9,
    response10$message(),AIC_resp10,response11$message(),AIC_resp11,
    response12$message(),AIC_resp12,file=filename,sep="\n",append=TRUE)

cat("\n",file=filename,append=TRUE)


cat("WITH OFFSET : \n",response7offset$message(),AIC_resp7offset,
    response8offset$message(),AIC_resp8offset,response9offset$message(),AIC_resp9offset,
    response10offset$message(),AIC_resp10offset,response11offset$message(),AIC_resp11offset,
    response12offset$message(),AIC_resp12offset,file=filename,sep="\n",append=TRUE)


